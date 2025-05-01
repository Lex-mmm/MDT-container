
from sdcProvider import SDCProvider
from datetime import datetime
import random
import time
client = SDCProvider()
client.start()

data = {'time': datetime.now(),  ## Round timestep to the second
        'values':{
            "mon_hr": 80,
            "mon_sat": random.randint(90, 100)
    } }

for i in range(100):
    data['time'] = datetime.now()
    data['values']['heart_rate'] = random.randint(60, 100)
    data['values']['SaO2'] = random.randint(90, 100)
    ## Send the data to the SDC
    client.sendData("vital_sign", data)
    time.sleep(1) ## Sleep for 1 second