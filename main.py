from flask import Flask, request, jsonify, render_template, redirect, url_for
import random
import json
import time
import calendar
import threading
import os
from digital_twin_model import DigitalTwinModel  # Import your updated class
from redis import Redis
from Comms.redisClient import RedisInit
import time, random
import calendar
import json

app = Flask(__name__, template_folder='templates', static_folder='static')

# Dictionary to store instances of DigitalTwinModel keyed by patient_id
twin_instances = {}

## Define global class for patient-level information
class Patient:

    def __init__(self, patient_id, param_file="parameters.json", pat_char=None, sleep=True):

        ## Initialize starting instance variables
        # Use provided database ID if available, otherwise generate random ID
        patient_db_id = os.environ.get('PATIENT_DB_ID')
        if patient_db_id and patient_db_id.strip():
            self.patient_id = patient_db_id.strip()
            print(f"📋 Using provided patient database ID: {self.patient_id}")
        else:
            self.patient_id = str(random.randint(1000, 9999))
            print(f"🎲 Generated random patient ID: {self.patient_id}")
        self.param_file = param_file
        self.running = False ## start right away
        self.data_epoch = []
        ## redis client for communication
        self.redisClient = RedisInit()
        # REMOVED: self.redisClient.redis.flushdb() - This was clearing ALL patients from Redis!
        # Only clear this specific patient's data if needed for restart
        try:
            # Clean up any existing data for THIS patient only
            existing_keys = self.redisClient.redis.keys(f"*{self.patient_id}*")
            if existing_keys:
                self.redisClient.redis.delete(*existing_keys)
                print(f"Cleaned up existing data for patient {self.patient_id}")
        except Exception as e:
            print(f"Note: Could not clean existing patient data: {e}")
        
        ## re-initialize the redis client
        self.redisClient = RedisInit()
        self.ingestDemographics()

        self.model = DigitalTwinModel(self.patient_id, param_file, data_callback=None, 
                                     sleep=sleep)
        
        self.model.redisClient = self.redisClient ## enable use further down the line
        
    def ingestDemographics(self):
        ## Ingest patient demographics into Redis
        # Use custom name if provided via environment, otherwise random
        patient_name = os.environ.get('PATIENT_NAME', '').strip()
        
        if patient_name and ' ' in patient_name:
            name_parts = patient_name.split(' ', 1)
            first_name = name_parts[0]
            last_name = name_parts[1]
        elif patient_name:
            first_name = patient_name
            last_name = random.choice(['Loon, van', 'Kappen', 'Koomen', 'Klein', 'Linden, van der'])
        else:
            first_name = random.choice(['Erik', 'Lex', 'Teus','Jan', 'Piet', 'Kees', 'Henk'])
            last_name = random.choice(['Loon, van', 'Kappen', 'Koomen', 'Klein', 'Linden, van der'])
        gender = random.choice(['M', 'F'])
        dob_year = random.randint(1950, 2000)
        dob_month = random.randint(1, 12)
        dob_day = random.randint(1, 28)
        dob = f"{dob_year}-{dob_month:02d}-{dob_day:02d}"
        dobNumeric = calendar.timegm(time.strptime(dob, "%Y-%m-%d")) ## convert to epoch time
        payload = {
            "ptID": str(self.patient_id),
            "first_name": str(first_name),
            "last_name": str(last_name),
            "dob": (dobNumeric),
            "gender": str(gender)
        }
        ## Ingest into Redis
        print(f"📝 Patient demographics: {payload}")
        patient_key = f"PATIENTS:{self.patient_id}"
        
        try:
            if not self.redisClient.redis.exists(patient_key):
                print(f"✨ Creating new patient {self.patient_id} in Redis")
                self.redisClient.redis.execute_command("JSON.SET", patient_key, "$", json.dumps(payload))
                print(f"✅ Patient {self.patient_id} created successfully in Redis")
            else:
                print(f"ℹ️  Patient {self.patient_id} already exists in Redis, updating...")
                self.redisClient.redis.execute_command("JSON.SET", patient_key, "$", json.dumps(payload))
                print(f"✅ Patient {self.patient_id} updated successfully in Redis")
                
            ## Verify the data was stored
            stored_data = self.redisClient.redis.execute_command("JSON.GET", patient_key)
            print(f"🔍 Verification - stored data: {stored_data}")
            
        except Exception as e:
            print(f"❌ Error storing patient demographics: {e}")
            
        try:
            ## Check the demographics index
            check = self.redisClient.redis.execute_command('FT.SEARCH', 'demographics_index', '*')
            print(f"📋 Demographics index check: {check}")
        except Exception as e:
            print(f"⚠️  Could not check demographics index: {e}")


    def eventListen(self):
        streamKey = f"EVENT:{self.patient_id}"
        print(f"🎧 Listening for events on stream: {streamKey}")
        
        # Track the last processed message ID to avoid re-processing
        last_id = '0'  # Start from beginning only on first read
        
        ## Create the stream if it doesn't exist
        if not self.redisClient.redis.exists(streamKey):
            self.redisClient.redis.xadd(streamKey, {'init': '1'}, id='*')
            print(f"✨ Stream {streamKey} created.")
        
        while self.running:
            try:
                # Read new messages from the last processed ID with blocking
                # Using block=1000 (1 second) to avoid busy waiting while still being responsive
                messages = self.redisClient.redis.xread({streamKey: last_id}, count=1, block=1000)
                
                if messages:
                    for stream, entries in messages:
                        for entry in entries:
                            message_id, message_data = entry  # Unpack the message ID and data
                            print(f"📨 Stream: {stream}, ID: {message_id}, Data: {message_data}")
                            
                            # Update last processed ID immediately to avoid re-processing
                            last_id = message_id
                            
                            # Process the message data
                            if message_data == {'init': '1'}:
                                print("🔧 Initialization message received. Ignoring.")
                                self.redisClient.redis.xdel(streamKey, message_id)  # Delete the message after processing
                                continue
                            else:
                                print(f"⚡ Processing therapy event: {message_data}")
                                
                                # Validate required fields before processing
                                if 'event' in message_data and 'eventSeverity' in message_data and 'eventType' in message_data:
                                    self.processEvent(
                                        message_data['event'], 
                                        message_data['eventSeverity'], 
                                        message_data['eventType']
                                    )
                                    print(f"✅ Successfully processed {message_data['event']} therapy for patient {self.patient_id}")
                                else:
                                    print(f"⚠️  Invalid event data format: {message_data}")
                                
                                # Delete the message after successful processing
                                self.redisClient.redis.xdel(streamKey, message_id)
                else:
                    # No new messages within the timeout period - this is normal
                    # Using shorter sleep since we're already using blocking read
                    time.sleep(0.1)
                    
            except Exception as e:
                print(f"❌ Error in event listening loop: {e}")
                # Sleep a bit before retrying to avoid tight error loops
                time.sleep(1)



    

    def startup(self):
        """Start the simulation."""
        self.running = True
        self.model.start_simulation()


    def stop_patient(self):
        """Stop the simulation."""
        self.running = False
        self.model.stop_simulation()

    def processEvent(self, event, eventSeverity, eventType=None):
        """Process event for patient consequences."""
        print(f"🔄 ----- PROCESSING EVENT ----- {event} ----- {eventSeverity} ----- {eventType}")
        
        if not event or eventSeverity is None:
            error_msg = f"❌ Missing data: event={event}, eventSeverity={eventSeverity}"
            print(error_msg)
            return error_msg
            
        processedEvent = None ## start with empty, fill if needed
        
        try:
            if eventType == "disease" or eventType == "recovery":
                print(f"🦠 Processing pathology event: {event}")
                processedEvent = self.model.pathologies.processPathology(event, eventSeverity)
            elif eventType == "treatment":
                print(f"💊 Processing therapy treatment: {event} with severity {eventSeverity}")
                processedEvent = self.model.therapeutic.processTherapeutic(event, eventSeverity)
            else:
                print(f"⚠️  Unknown event type: {eventType}")
                return f"Unknown event type: {eventType}"
            
            print(f"🧬 Processed event result: {processedEvent}")
            
            if processedEvent:
                for processed_item in processedEvent:
                    ## Allow multiple actions to flow from a single event (e.g. starting point + decay)
                    self.model.addProcessedEvent(processed_item)
                    print(f"✅ Added processed event to model: {processed_item.get('event', 'unknown')}")
            else:
                print(f"⚠️  No processed event returned for {event} (this might be expected for some events)")
                
        except Exception as e:
            error_msg = f"❌ Error processing event {event}: {e}"
            print(error_msg)
            return error_msg

        ## Single-event-layout:
        # event: {
        #     "event": "disease",
        #     "eventSeverity": 2,
        #     "eventType": "disease"
        #     "timeCategorical": "continuous" || "limited"
        #     "lastEmission": self.dt
        #     "timeInterval": 0.5
        #     "timeUnit": "hours"
        #     "eventCount": 5
        #     "parameters": {name: x, value: y, action: 'set' vs 'relative', type: 'absolute' vs 'relative}
        #}
            




### START THE SERVER ###
"""Start simulation for a specific patient."""
## Get patient configuration from environment variables
# Patient ID will be randomly generated by the Patient class
container_id = os.environ.get('PATIENT_ID', f"Container-{random.randint(1000, 9999)}")  # Used for container tracking
param_file = os.environ.get('PARAM_FILE', 'healthyFlat.json')
sleep = os.environ.get('SLEEP', 'false').lower() == 'true'

print(f"🏥 Starting patient container:")
print(f"   Container ID: {container_id}")
print(f"   Parameter File: {param_file}")
print(f"   Patient Name: {os.environ.get('PATIENT_NAME', 'Auto-generated')}")
print(f"   Sleep Mode: {sleep}")

if container_id not in twin_instances:
    # Create a new patient instance (patient ID will be randomly generated)
    twin_instances[container_id] = Patient(patient_id= None, param_file= param_file, sleep=sleep)


patient = twin_instances[container_id]
if not patient.running:
    # Start the simulation in a new thread
    print(f"🚀 Starting simulation threads for patient {patient.patient_id} (container: {container_id})")
    threading.Thread(target=patient.startup).start()
    threading.Thread(target=patient.eventListen).start()
    print(f"✅ Simulation started successfully for patient {patient.patient_id}")
    print(f"🏥 Patient container is now running and ready for monitoring!")
else:
    print(f"ℹ️  Simulation already running for patient {patient.patient_id}")
