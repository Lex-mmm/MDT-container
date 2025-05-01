
import json


class Therapy:
    """
    Base class for all therapies.
    """

    def __init__(self):
        with open("healthyFlat.json", "r") as f:
            self.parameters = json.load(f)

        with open("Events/Therapy/therapy.json", "r") as f:
            self.therapy = json.load(f)
    

    def processTherapeutic(self, event, eventSeverity):
        """
        Process the therapeutic event and update the master parameters.
        """

        therapy = self.therapy[event]
        if therapy:
            therapyData = therapy[f"severity_{eventSeverity}"]['intervention']['parameters']
            for paramName in therapyData:
                ## Name of the parameter to alter
                value = therapyData[paramName]['value']
                type = therapyData[paramName]['type']
                rate = therapyData[paramName]['rate']
                action = therapyData[paramName]['action']
                
                rateTimeValue = rate['timeValue']
                rateTimeUnit = rate['timeUnit']
                rateTimeCount = rate['count']
                timeCategorical = 'continuous' if isinstance(rateTimeCount, str) else 'limited'
                startShell = {
                    "event": event,
                    "eventSeverity": eventSeverity,
                    "eventType": "common",
                    "timeCategorical": timeCategorical,
                    "lastEmission": 0,
                    "timeInterval": rateTimeValue,
                    "timeUnit": rateTimeUnit,
                    "eventCount": rateTimeCount,  # Initial use only
                    "parameters": []
                }
                paramFillShell = {
                    "name": paramName,
                    "value": value,
                    "action": action,
                    "type": type
                }
                startShell["parameters"].append(paramFillShell)
            
            print(f"Starting shell is: {startShell}")
            return [ startShell ]
            