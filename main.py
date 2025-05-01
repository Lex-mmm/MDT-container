from flask import Flask, request, jsonify, render_template, redirect, url_for
import random
import json
import time
import calendar
import threading
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
        self.patient_id = random.randint(1000, 9999) ## random patient ID
        self.param_file = param_file
        self.running = False ## start right away
        self.data_epoch = []
        ## redis client for communication
        self.redisClient = RedisInit()
        self.redisClient.redis.flushdb() ## clear the database for testing purposes
        ## re-initialize the redis client
        self.redisClient = RedisInit()
        self.ingestDemographics()

        self.model = DigitalTwinModel(self.patient_id, param_file, data_callback=None, 
                                     sleep=sleep)
        
        self.model.redisClient = self.redisClient ## enable use further down the line
        
    def ingestDemographics(self):
        ## Ingest patient demographics into Redis
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
        print(f"Patient demographics: {payload}")
        if not self.redisClient.redis.exists(f"PATIENTS:{self.patient_id}"):
            #print(f"Creating patient {self.patient_id} in Redis.")
            self.redisClient.redis.execute_command("JSON.SET", f'PATIENTS:{self.patient_id}', "$", json.dumps(payload))
            ## Create the index for the patient demographics
            check = self.redisClient.redis.execute_command(
                    'FT.SEARCH', 'demographics_index', '*'
                )
            print(f"Patient demographics created: {check}")


    def eventListen(self):
        streamKey = f"EVENT:{self.patient_id}"
        print(f"Listening for events on stream: {streamKey}")
        ## Create the stream if it doesn't exist
        if not self.redisClient.redis.exists(streamKey):
            self.redisClient.redis.xadd(streamKey, {'init': '1'}, id='*')
            print(f"Stream {streamKey} created.")
        while self.running:
            messages = self.redisClient.redis.xread({streamKey: '0'}, count=1, block=0)
            if messages:
                for stream, entries in messages:
                    for entry in entries:
                        message_id, message_data = entry  # Unpack the message ID and data
                        print(f"Stream: {stream}, ID: {message_id}, Data: {message_data}")
                        
                        # Process the message data as needed
                        # Example: Convert JSON string to dictionary
                        
                        if message_data == {'init': '1'}:
                            print("Initialization message received. Ignoring.")
                            self.redisClient.redis.xdel(streamKey, message_id)  # Delete the message after processing

                            continue
                        else:
                            print(f"Processing: {message_data}")
                            self.processEvent(message_data['event'], message_data['eventSeverity'], message_data['eventType'])
                            self.redisClient.redis.xdel(streamKey, message_id)  # Delete the message after processing
            else:
                print("No new messages. Sleeping for a bit.")
                time.sleep(5)



    

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
        print(f" ----- PROCESSING EVENT ----- {event} ----- {eventSeverity} ----- {eventType}")
        if not event or not eventSeverity:
            return "Missing data: event or eventSeverity"
        processedEvent = None ## start with empty, fill if needed
        if eventType == "disease" or eventType == "recovery":
            processedEvent = self.model.pathologies.processPathology(event, eventSeverity)
        elif eventType == "treatment":
            processedEvent = self.model.therapeutic.processTherapeutic(event, eventSeverity)
        
        print(f"Processed event: {processedEvent}")
        if processedEvent:
            for event in processedEvent:
                ## Allow multiple actions to flow from a single event (e.g. starting point + decay)
                self.model.addProcessedEvent(event)

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
## Allow for back-end submission of patient_id and param_file
patient_id = "Erik"
param_file = "healthyFlat.json"
sleep = False

if patient_id not in twin_instances:
    # Create a new patient instance
    twin_instances[patient_id] = Patient(patient_id= patient_id, param_file= param_file, sleep=sleep)


patient = twin_instances[patient_id]
if not patient.running:
    # Start the simulation in a new thread
    threading.Thread(target=patient.startup).start()
    threading.Thread(target=patient.eventListen).start()
    print(f"Simulation started for patient {patient_id}")
#    return jsonify({"status": f"Simulation started for patient {patient_id}"})
#else:
 #   return jsonify({"status": f"Simulation already running for patient {patient_id}"})
