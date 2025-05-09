from PyQt6 import QtWidgets  # Should work with PyQt5 / PySide2 / PySide6 as well
import pyqtgraph as pg
from PyQt6.QtWidgets import QApplication, QMainWindow, QLabel
from PyQt6.QtCore import QTimer, Qt
from random import sample
import requests
import pyodbc
from datetime import datetime, timedelta
# Always start by initializing Qt (only once per application)

from redis import Redis
import docker
import json
import time

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        

        self.redis = Redis(host='localhost', port=6379, db=0, decode_responses=True)
        self.redis.ping()

               
        # Define a top-level widget to hold everything
        self.w = QtWidgets.QWidget()

        # Create some widgets to be placed inside
        self.alarm_list = QtWidgets.QListWidget(  )
        self.pt_select = QtWidgets.QComboBox()
        self.pt_select.currentTextChanged.connect(self.update_selection) ## Refresh dependent on active patients

        self.plot_pt = pg.plot(
            title='Patient data',
            labels={'left': ('Value'), 'bottom': ('Time (s) ')} )
        self.plot_pt.addLegend()
        self.plot_pt.setYRange(0, 130)

        self.firstEvent = QtWidgets.QPushButton('MI-1')
        self.firstEvent.clicked.connect(self.event_handler)

        self.secEvent = QtWidgets.QPushButton('MI-2')
        self.secEvent.clicked.connect(self.event_handler)

        self.thirdEvent = QtWidgets.QPushButton('MI-3')
        self.thirdEvent.clicked.connect(self.event_handler)

        self.fourthEvent = QtWidgets.QPushButton('FI_O2 | 0.1')
        self.fourthEvent.clicked.connect(self.event_handler)

        self.fifthEvent = QtWidgets.QPushButton('O2 | 5l/min')
        self.fifthEvent.clicked.connect(self.event_handler)

        self.sixthEvent = QtWidgets.QPushButton('NaCl | 500CC')
        self.sixthEvent.clicked.connect(self.event_handler)



        

        # Create a grid layout to manage the widgets size and position
        self.layout = QtWidgets.QGridLayout()

        # Add widgets to the layout in their proper positions
        self.layout.addWidget(self.alarm_list, 2, 0)  # list widget goes in bottom-left
        self.layout.addWidget(self.plot_pt, 0, 1, 5, 1)  # plot goes on right side, spanning 3 rows
        self.layout.addWidget(self.pt_select, 0, 0)
        ## add another row to the far left with multiple buttons
        self.layout.addWidget(self.firstEvent, 0, 7)  # combo box goes in top-left
        self.layout.addWidget(self.secEvent, 1, 7)
        self.layout.addWidget(self.thirdEvent, 2, 7)
        self.layout.addWidget(self.fourthEvent, 3, 7)
        self.layout.addWidget(self.fifthEvent, 4, 7)
        self.layout.addWidget(self.sixthEvent, 5, 7)
        


        self.w.setLayout(self.layout)
        self.setCentralWidget(self.w)

        self.refresh_plots = QTimer(self)
        self.refresh_plots.timeout.connect(self.update_plot)
        self.refresh_plots.start(1000)

        self.refresh_data = QTimer(self)
        self.refresh_data.timeout.connect(self.update_data)
        self.refresh_data.start(1000)

        self.refresh_patients = QTimer(self)
        self.refresh_patients.timeout.connect(self.update_patients)
        self.refresh_patients.start(2500)

        self.refresh_alarms = QTimer(self)
        self.refresh_alarms.timeout.connect(self.update_alarms)
        self.refresh_alarms.start(1000)

        self.data_timestamp = []
        self.data_mon_sat = []
        self.data_mon_hr = []
        self.data_mon_map = []
        self.data_mon_rr = []

        self.pat_id = ''

    def event_handler(self):
            sender = self.sender().text()
            match sender:
                case 'MI-1':
                    payload = {
                        'event': 'myocardialInfarction',
                        'eventSeverity': '1',
                        'eventType': 'disease'
                    }
                    print("Event triggered: MI-1")
                case 'MI-2':
                    payload = {
                        'event': 'myocardialInfarction',
                        'eventSeverity': '2',
                        'eventType': 'disease'
                    }
                    print("Event triggered: MI-2")

                case 'MI-3':
                    payload = {
                        'event': 'myocardialInfarction',
                        'eventSeverity': '3',
                        'eventType': 'disease'
                    }
                    print("Event triggered: MI-3")
                case 'O2 | 5l/min':
                    payload = {
                        'event': 'O2-Suppletion',
                        'eventSeverity': '5',
                        'eventType': 'treatment'
                    }
                    print("Event triggered: O2 | 5l/min")

                case 'NaCl | 500CC':
                    payload = {
                        'event': 'NaCl-suppletion',
                        'eventSeverity': '500',
                        'eventType': 'treatment'
                    }
                    print("Event triggered: NaCl suppletion")
                case 'FI_O2 | 0.1':
                    payload = {
                        'event': 'O2-Suppletion',
                        'eventSeverity': '0',
                        'eventType': 'treatment'
                    }
                    print("Event triggered: FI_O2 | 0.1")
                case _:
                    print(f"Unknown event: {sender}")
            # Send the payload to the server
            stream_key = f"EVENT:{self.pat_id}"
            self.redis.xadd(stream_key, payload, id='*')



    def update_selection(self, choice):
        choiceNumeric = choice.split('   |   ')[0]  # Extract the ptID from the selected item
        self.pat_id = choiceNumeric

    def update_alarms(self):
            if self.pat_id == '' or self.pat_id == None:
                pass
            else:

                result = self.redis.execute_command(
                    'FT.SEARCH', 'alarm_index',
                    f'@ptID:{{{str(self.pat_id)}}}',  # Query by ptID
                    'SORTBY', 'timestamp', 'ASC',
                    'LIMIT', '0', '25'  # Get the last 10 entries
                )

                ## Add the alarms to the list if not already present
                if result and len(result) > 1:
                    data = self.deserialize(result)
                    #print(f"Alarm for patient {self.pat_id}: {result}")
                    #return result
                    if data and len(data) > 1:
                        #print(f"Alarm for patient {self.pat_id}: {data}")
                        for entry in data:
                            entryData = data[entry]
                            #print(entryData)
                            timeFormatted = datetime.fromtimestamp(float(entryData['timestamp'])).strftime('%H:%M:%S')
                            alarm = f"{timeFormatted} \t  {entryData['alarm_status']} \t {entryData['alarm_msg']}"

                            if alarm not in [self.alarm_list.item(i).text() for i in range(self.alarm_list.count())]:
                                self.alarm_list.insertItem(0, alarm)
            
    def update_patients(self):
        result = self.redis.execute_command(
            'FT.AGGREGATE', 'demographics_index',
            '*',  # Match all entries
            'GROUPBY', '1', '@ptID'
        )

        # Parse the result to extract unique ptIDs
        unique_ptIDs = [row[1] for row in result[1:]]
        for selectID in unique_ptIDs:
            ID_data = self.redis.execute_command(
                'FT.SEARCH', 'demographics_index',
                f'@ptID:{{{str(selectID)}}}'
                )
            ID_data = self.deserialize(ID_data)[f'PATIENTS:{selectID}']
            ID_line = f'{ID_data["ptID"]}   |   {ID_data["last_name"]},  {ID_data["first_name"]}'
            current_items = [self.pt_select.itemText(i) for i in range(self.pt_select.count())]
            if ID_line not in current_items:
                self.pt_select.addItem(ID_line)

    def deserialize(self, data):
        ## input is key-value of options as resulted from the search-index
        ## output is a list of dictionaries with the key-value pairs
        deserialized = {}
        if not data or len(data) < 2:
            print("No data found or invalid data format.")
            return None
        data = data[1:]  # Skip the first element which is the count
        
        for i in range(0, len(data), 2):

            # key = name of the input, can be ignored -> Future = reference number unique to that entry
            key = data[i]
            
            # FUTURE: Implement below
            #json_data = results[i+1][ int([k for k in range(len(results[i+1])) if results[i+1][k] == '$' ] + 1 )  ] ## '$'-sign +1 = actual data, before is sorting data and index-value 
            #print(data[i+1])
            json_data = data[i+1][len(data[i+1]) -1 ] ## last element is the actual data
            #print(json_data)
            #data = json.loads(json_data)
            #data['DATUM'] = time.strftime('%Y-%m-%d', time.gmtime(data['DATUM']))
            deserialized[key] = json.loads(json_data)

        return deserialized
    
    def update_data(self):
        if self.pat_id == '' or self.pat_id is None:
            return
        
        ### QUERYING PREVIOUS VITALS
        try:
            data = self.redis.execute_command(
                'FT.SEARCH', 'vital_sign_index',
                f'@ptID:{{{self.pat_id}}}',  # Query by ptID
                'SORTBY', 'timestamp', 'DESC',
                'LIMIT', '0', '120'  # Get the last 10 entries
            )
            #print(f"Vital signs for patient {self.pat_id}: {result}")
            #return result
        except Exception as e:
            print(f"Error fetching vital signs: {e}")
            #return []
        
        if data:
            data = self.deserialize(data)

            timestamp = []
            mon_sat = []
            mon_hr = []
            mon_map = []
            mon_rr = []

            for entry in data:
                #try:
                    #print(f"Data type: {type(entry['values']['time'])}")
                    entryData = data[entry]

                    timestamp.append(datetime.fromtimestamp((float(entryData['timestamp']))))
                    #timestamp.append(datetime.strptime(float(entryData['timestamp']), '%a, %d %b %Y %H:%M:%S %Z'))
                    mon_sat.append(float(entryData['mon_sat']))
                    mon_hr.append(float(entryData['mon_hr']))
                    mon_map.append(float(entryData['mon_ibp_mean']))
                    mon_rr.append(float(entryData['mon_rr']))
                #except (ValueError, TypeError, KeyError):
                    # Skip invalid entries
                    #print(f"Invalid entry: {entry}")
                    #continue

            self.data_timestamp = [(time - max(timestamp)).total_seconds() for time in timestamp]  # Normalize timestamps to lead up to t=0
            self.data_mon_sat = mon_sat
            self.data_mon_hr = mon_hr
            self.data_mon_map = mon_map
            self.data_mon_rr = mon_rr
            

    def update_plot(self):

        if self.pat_id == '' or self.pat_id == None:
            pass
        else:
            self.plot_pt.clear()
            self.plot_pt.plot(name ="SaO2", x = self.data_timestamp, y = self.data_mon_sat, pen = pg.mkPen(color = 'b', width = 3))
            self.plot_pt.plot(name = "HR", x = self.data_timestamp, y = self.data_mon_hr, pen = pg.mkPen(color = 'g', width = 3))
            self.plot_pt.plot(name = "MAP", x = self.data_timestamp, y = self.data_mon_map, pen = pg.mkPen(color = 'r', width = 3))
            self.plot_pt.plot(name = "RR", x = self.data_timestamp, y = self.data_mon_rr, pen = pg.mkPen(color = 'y', width = 3))
            self.plot_pt.addLegend()


if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    window = MainWindow()
    window.show()
    app.exec()

