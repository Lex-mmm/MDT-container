
from redis import Redis
import time
#import docker
import json
import calendar

class RedisInit:
    def __init__(self, host='host.docker.internal', port=6379):
        ## startup the docker if not yet started

        ## Only necessary if running locally
        '''
        
        # Initialize Docker client
        try:
            client = docker.from_env()
            print("Docker client initialized")
        except docker.errors.DockerException as e:
            print(f"Error initializing Docker client: {e}")
            exit(1)

        container_name = 'redis-stack'
        try:
            container = client.containers.get(container_name)
            if container.status != 'running':
                container.start()
                time.sleep(2)
            else: 
                container.kill()
                time.sleep(2)
                container.start()
                time.sleep(2)
            print("RedisJSON Docker container is running")
        except docker.errors.NotFound:
            print("RedisJSON Docker container not found. Creating and starting a new one...")
            container = client.containers.run('redis/redis-stack:latest', name=container_name,
                                            ports={'6379/tcp': 6379}, detach=True)
            time.sleep(2)
            print("RedisJSON Docker container started")

        '''
        

        self.redis = Redis(host=host, port=port,  decode_responses=True)
        #self.redis.ping()  # Check if the connection is successful
        print("Connected to Redis")
        self.createIndices()
        print("Indices created")
 

    def createIndices(self):

        
        try:
            print("Initializing vital signs index...")
            self.redis.execute_command(
                    'FT.CREATE', 'vital_sign_index', 'ON', 'JSON', 'PREFIX', '1', 'VITALS:',
                    'SCHEMA',
                    '$.ptID', 'AS', 'ptID', 'TAG',
                    '$.timestamp', 'AS', 'timestamp', 'NUMERIC', 'SORTABLE',
                    '$.mon_hr', 'AS', 'mon_hr', 'NUMERIC', 'SORTABLE',
                    '$.mon_sat', 'AS', 'mon_sat', 'NUMERIC', 'SORTABLE',
                    '$.mon_ibp_dia', 'AS', 'mon_ibp_dia', 'NUMERIC', 'SORTABLE',
                    '$.mon_ibp_sys', 'AS', 'mon_ibp_sys', 'NUMERIC', 'SORTABLE',
                    '$.mon_ibp_mean', 'AS', 'mon_ibp_mean','NUMERIC', 'SORTABLE',
                    '$.mon_rr', 'AS', 'mon_rr','NUMERIC', 'SORTABLE',
            )
            
            print("Vital signs index created.")
        except Exception as e:
            print(f"Error creating indices: {e}")
            print("VITAL SIGN INDEX.")
        
        try:
            print("Initializing alarm index...")
            self.redis.execute_command(
                    'FT.CREATE', 'alarm_index', 'ON', 'JSON', 'PREFIX', '1', 'ALARMS:',
                    'SCHEMA',
                    '$.ptID', 'AS', 'ptID', 'TAG',
                    '$.timestamp', 'AS', 'timestamp', 'NUMERIC', 'SORTABLE',
                    '$.alarm_msg', 'AS', 'alarm_msg', 'TEXT',
                    '$.alarm_prio', 'AS', 'alarm_prio', 'TAG',
                    '$.alarm_status', 'AS', 'alarm_status', 'TAG',
            )
            print("Alarm index created.")
        except Exception as e:
            print(f"Error creating indices: {e}")
            print("ALARM INDEX.")

        ## Patient demographics
        try:
            print("Initializing patient index...")
            self.redis.execute_command(
                    'FT.CREATE', 'demographics_index', 'ON', 'JSON', 'PREFIX', '1', 'PATIENTS:',
                    'SCHEMA',
                    '$.ptID', 'AS', 'ptID', 'TAG',
                    '$.first_name', 'AS', 'first_name', 'TEXT',
                    '$.last_name', 'AS', 'last_name', 'TEXT',
                    '$.dob', 'AS', 'dob', 'NUMERIC', 'SORTABLE',
                    '$.gender', 'AS', 'gender', 'TAG',
            )
            print("patient index created.")
        except Exception as e:
            print(f"Error creating indices: {e}")
            print("DEMOGRAPHICS INDEX.")


    def add_vital_sign(self, ptID, payload):
        #print(f"Adding vital sign for patient {ptID} with payload: {payload}")
        vitalEntry = {
            "ptID": ptID,
            "timestamp": int(calendar.timegm(payload['time'].timetuple())),
            "mon_hr": float(payload['values']["heart_rate"]),
            "mon_sat": float(payload['values']["SaO2"]),
            "mon_ibp_dia": float(payload['values']["DAP"]),
            "mon_ibp_sys": float(payload['values']["SAP"]),
            "mon_ibp_mean": float(payload['values']["MAP"]),
            "mon_rr": float(payload['values']["RR"]),
        }
        #print(f"Vital sign entry: {vitalEntry}")
        ## insert the entry into the VITALS index table:
        try:
            key = f"VITALS:{ptID}:{int(vitalEntry['timestamp'])}"  # Ensure key matches the prefix
            self.redis.execute_command("JSON.SET", key, "$", json.dumps(vitalEntry))
            self.sdcComms("VITALS", payload) ## Send over SDC

        except Exception as e:
            print(f"Error inserting vital sign entry: {e}")

            
    def add_alarm(self, ptID, payload):
        payload['timestamp'] = int(calendar.timegm(payload['timestamp'].timetuple()))
        print(f"Adding alarm for patient {ptID} with payload: {payload}")
        ## insert the entry into the ALARMS index table:
        try:
            key = f"ALARMS:{payload['ptID']}:{payload['timestamp']}"  # Ensure key matches the prefix
            result = self.redis.execute_command("JSON.SET", key, "$", json.dumps(payload))
            print(f"Result: {result}")
            self.sdcComms("ALARMS", payload) ## Send over SDC
            

        except Exception as e:
            print(f"Error adding alarm entry: {e}")
    


    def sdcComms(self, eventType, payload):
        ## send the data to the SDC
            # EventType = "VITALS" or "ALARMS"

        pass