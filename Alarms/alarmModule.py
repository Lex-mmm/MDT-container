import json
from Alarms.AscomMathFunctions import AscomMathFunctions
from datetime import datetime, timedelta
from Comms.redisClient import RedisInit

class alarmModule:
    def __init__(self, pt_id):
        with open ("Alarms/static/initialSetting.json", "r") as f:
            payload = json.load(f)
            self.thresholdAlarm = payload["AlarmSettingThreshold"]
            self.algoAlarm = payload["AlarmSettingAlgo"]
            self.algoTightness = payload["AlgoTightness"]
            self.algoReferenceRange = payload["ReferenceRangeValues"]
        
        with open ("Alarms/static/messageConversion.json", "r") as f:
            payload = json.load(f)
            self.messageConversion = payload["AlarmMessage"]
        
        self.AlgoFunctions = AscomMathFunctions()
        self.AlgoHorizon = 60 ## 60-seconds horizon for the algorithm-based alarms

        self.alarms = []
        self.comms = RedisInit()
        self.ptID = pt_id

    def evaluate_data(self, curr_data, historic_data):
        ''' Evaluate the alarm state for the current data.
            # Format:
            # { data = {'time': np.round(self.t),  ## Round timestep to the second
                        'values':{
                            "heart_rate": np.round(self.current_heart_rate, 2),
                            "SaO2": np.round(self.current_SaO2, 2),
                            "MAP": np.round(np.mean(self.P_store), 2),
                            "RR": np.round(self.RR, 2),
                            "etCO2": np.round(self.current_state[17], 2)
                    } }
                    }
                    '''
        ## Create timeframe-specific set based on Algo-horizon
        curr_time = curr_data["time"]
        start_time = curr_time - timedelta(seconds = self.AlgoHorizon)
        end_time = curr_time
        AlgoDataTime = []
        AlgoDataValSaO2 = []
        AlgoDataValHR = []
        AlgoDataValMAP = []
        AlgoDataValRR = []


        for data in historic_data:
            if start_time <= data["time"] <= end_time:
                AlgoDataTime.append(data["time"])
                AlgoDataValSaO2.append(data["values"]["SaO2"])
                AlgoDataValHR.append(data["values"]["heart_rate"])
                AlgoDataValMAP.append(data["values"]["MAP"])
                AlgoDataValRR.append(data["values"]["RR"])

        ## Oxygen saturation alarms:
        self.ALGO_checkAlarmValue({'time': AlgoDataTime, 'values':AlgoDataValSaO2}, "OxygenSaturation")
        self.THRESHOLD_checkAlarmValue(curr_data["time"], curr_data["values"]["SaO2"], "OxygenSaturation")
        ## Heart rate alarms:
        self.ALGO_checkAlarmValue({'time': AlgoDataTime, 'values':AlgoDataValHR}, "HeartRate")
        self.THRESHOLD_checkAlarmValue(curr_data["time"], curr_data["values"]["heart_rate"], "HeartRate")
        ## Mean arterial pressure alarms:
        self.ALGO_checkAlarmValue({'time': AlgoDataTime, 'values':AlgoDataValMAP}, "BloodPressureMean")
        self.THRESHOLD_checkAlarmValue(curr_data["time"], curr_data["values"]["MAP"], "BloodPressureMean")
        ## Respiratory rate alarms:
        self.ALGO_checkAlarmValue({'time': AlgoDataTime, 'values':AlgoDataValRR}, "RespiratoryRate")
        self.THRESHOLD_checkAlarmValue(curr_data["time"], curr_data["values"]["RR"], "RespiratoryRate")
        



## Old, threshold-based system
    def getThresholdSetting(self, parameter):
        return self.thresholdAlarm[parameter]
    
    def evaluateThresholdAlarmState(self, alarmType, parameter, state, message, timestamp, priority):
        curr_state = self.thresholdAlarm[parameter]["AlarmState"][alarmType]
        if curr_state != state and state == True: # alarm state changed, activated state
            self.thresholdAlarm[parameter]["AlarmState"][alarmType] = state
            self.alarmTrigger(parameter, alarmType, message, timestamp, priority)

        elif curr_state != state and state == False: # alarm state changed, deactivated
            self.thresholdAlarm[parameter]["AlarmState"][alarmType] = state
            self.alarmResolve(parameter, alarmType, message, timestamp, priority)
    

    def alarmTrigger(self, parameter, alarmType, message, timestamp, priority):  
        ## trigger callback function
        payload = {
            "ptID": str(self.ptID),
            "timestamp": timestamp,
            "alarm_msg": message,
            "alarm_prio": priority,
            "alarm_status": "ON",
            "parameter": parameter,
            "alarm_type": alarmType
        }
        self.comms.add_alarm(self.ptID, payload)

        self.alarms.insert(0, f"{payload['timestamp']} | TRIGGERED: {payload['alarm_msg']}")
        #print(f"Alarm TRIGGERED for {parameter} with type {alarmType}. Message: {message}")

    def alarmResolve(self, parameter, alarmType, message, timestamp, priority):
        ## resolve callback function
        payload = {
            "ptID": str(self.ptID),
            "timestamp": timestamp,
            "alarm_msg": message,
            "alarm_prio": priority,
            "alarm_status": "OFF",
            "parameter": parameter,
            "alarm_type": alarmType
        }
        self.comms.add_alarm(self.ptID, payload)
        self.alarms.insert(0, f"{payload['timestamp']} | resolved: {payload['alarm_msg']}")
        #print(f"Alarm RESOLVED for {parameter} with type {alarmType}.")

    def THRESHOLD_checkAlarmValue(self, timestamp, value, parameter):
        paramSetting = self.getThresholdSetting(parameter)
        if paramSetting["Active"] == True:
            ## pre-compile messages and priorities
            messageHigh, priorityHigh = self.matchAlarmMessage("HIGH", parameter)
            messageLow, priorityLow = self.matchAlarmMessage("LOW", parameter)

            if paramSetting["LowerLimit"] < value < paramSetting["UpperLimit"]:
                self.evaluateThresholdAlarmState("HIGH", parameter, False, messageHigh, timestamp, None)
                self.evaluateThresholdAlarmState("LOW", parameter, False, messageLow, timestamp, None)
                return None ## no alarm, all within range
            
            elif value > paramSetting["UpperLimit"]:
                self.evaluateThresholdAlarmState("HIGH", parameter, True, messageHigh, timestamp, priorityHigh)
            elif value < paramSetting["LowerLimit"]:
                self.evaluateThresholdAlarmState("LOW", parameter, True, messageLow, timestamp, priorityLow)
        
    def matchAlarmMessage(self, alarmType, parameter):
        if parameter not in self.messageConversion:
            return "UNKNOWN PARAMETER"
        else:
            message = self.messageConversion[parameter][alarmType]
            return message["message"], message["messagePriority"] if message else f"No message found for {alarmType} on {parameter}"
       

## New, algorithm-based system
    def getAlgoSetting(self, parameter):
        return self.algoAlarm[parameter]
    
    def getAlgoTightnessValue(self, tightness):
        setting = self.algoTightness[tightness]
        return setting["yellow"], setting["orange"], setting["red"]

    def ALGO_checkAlarmValue(self, data, parameter):
        paramSetting = self.getAlgoSetting(parameter)

        if paramSetting["Active"] == True:
            tightness = paramSetting["tightness"]
            yellow, orange, red = self.getAlgoTightnessValue(tightness)
            envelopEnsemble = self.AlgoFunctions.MATH_envelop(data, 
                                                           referenceValMax = self.algoReferenceRange[parameter]["UpperLimit"], 
                                                           referenceValMin = self.algoReferenceRange[parameter]["LowerLimit"])
            if envelopEnsemble:
                envelopValue = envelopEnsemble["envelop_total"] if envelopEnsemble["envelop_total"] else 0
            else:
                envelopValue = 0.0 ## no data outside of range
            
            ## pre-compile messages and priorities
            messageYellow, priorityYellow = self.matchAlarmMessage("ENVELOP_YELLOW", parameter)
            messageOrange, priorityOrange = self.matchAlarmMessage("ENVELOP_ORANGE", parameter)
            messageRed, priorityRed = self.matchAlarmMessage("ENVELOP_RED", parameter)

            if envelopValue >= red:
                self.evaluateAlgoAlarmState("red", parameter, True, messageRed, max(data["time"]), priorityRed)

            elif envelopValue >= orange:
                self.evaluateAlgoAlarmState("red", parameter, False, messageRed, max(data["time"]), None)
                self.evaluateAlgoAlarmState("yellow", parameter, False, messageYellow, max(data["time"]), None)
                ## orange alarm is active
                self.evaluateAlgoAlarmState("orange", parameter, True, messageOrange, max(data["time"]), priorityOrange)
            elif envelopValue >= yellow:
                self.evaluateAlgoAlarmState("red", parameter, False, messageRed, max(data["time"]), None)
                self.evaluateAlgoAlarmState("orange", parameter, False, messageOrange, max(data["time"]), None)
                self.evaluateAlgoAlarmState("yellow", parameter, True, messageYellow, max(data["time"]), priorityYellow)

            else: ## envelopValue < yellow - deactivate all alarms
                self.evaluateAlgoAlarmState("yellow", parameter, False, messageYellow, max(data["time"]), None)
                self.evaluateAlgoAlarmState("orange", parameter, False, messageOrange, max(data["time"]), None)
                self.evaluateAlgoAlarmState("red", parameter, False, messageRed, max(data["time"]), None)
                

    def evaluateAlgoAlarmState(self, alarmType, parameter, state, message, timestamp, priority):
        curr_state = self.algoAlarm[parameter]["AlarmState"][alarmType]
        if curr_state != state and state == True: # alarm state changed, activated state
            self.algoAlarm[parameter]["AlarmState"][alarmType] = state
            self.alarmTrigger(parameter, alarmType, message, timestamp, priority)

        elif curr_state != state and state == False: # alarm state changed, deactivated
            self.algoAlarm[parameter]["AlarmState"][alarmType] = state
            self.alarmResolve(parameter, alarmType, message, timestamp, priority)