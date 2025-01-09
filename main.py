from flask import Flask, request, jsonify, render_template, redirect, url_for
from flask_socketio import SocketIO, emit, join_room
import threading
from digital_twin_model import DigitalTwinModel  # Import your updated class
from Alarms.alarm_check import Alarms  # Import alarm-generation class
import json

app = Flask(__name__, template_folder='templates', static_folder='static')
socketio = SocketIO(app)

# Dictionary to store instances of DigitalTwinModel keyed by patient_id
twin_instances = {}

## Define global class for patient-level information
class Patient:

    def __init__(self, patient_id, param_file="parameters.json", pat_char=None):

        ## Initialize starting instance variables
        self.patient_id = patient_id
        self.param_file = param_file
        self.running = False
        self.data_epoch = []
        self.model = DigitalTwinModel(patient_id, param_file, data_callback=None, 
                                    alarm_callback = Alarms(self.patient_id))
        self.alarms = []

    def startup(self):
        """Start the simulation."""
        self.running = True
        self.model.start_simulation()


    def stop_patient(self):
        """Stop the simulation."""
        self.running = False
        self.model.stop_simulation()




@app.route('/start_simulation', methods=['POST'])
def start_simulation():
    """Start simulation for a specific patient."""
    ## Allow for back-end submission of patient_id and param_file
    patient_id = request.args.get("patient_id")
    param_file = request.args.get("param_file")

    if not patient_id:
        return jsonify({"error": "Missing patient ID"}), 400

    if patient_id not in twin_instances:
        # Create a new patient instance
        twin_instances[patient_id] = Patient(patient_id, param_file)
    

    patient = twin_instances[patient_id]
    if not patient.running:
        # Start the simulation in a new thread
        threading.Thread(target=patient.startup).start()
        return jsonify({"status": f"Simulation started for patient {patient_id}"})
    else:
        return jsonify({"status": f"Simulation already running for patient {patient_id}"})
    

@app.route('/stop_simulation', methods=['POST'])
def stop_simulation():
    """Stop simulation for a specific patient."""
    patient_id = request.args.get("patient_id")
    if patient_id in twin_instances and twin_instances[patient_id].running:
        twin_instances[patient_id].stop_patient()
        return jsonify({"status": f"Simulation stopped for patient {patient_id}"})
    else:
        return jsonify({"status": f"No running simulation for patient {patient_id}"})

## Call for getting currently active patients
@app.get(('/get_patient_list'))
def get_patient_list():
    """Get the list of patients."""
    return jsonify(list(twin_instances.keys()))

## Call for changing parameter settings -> Ter check over methodiek
@app.route('/set_param', methods=['POST'])
def set_param():
    """Change a parameter value for a specific patient."""
    patient_id = request.args.get("patient_id")
    param = request.args.get("param")
    value = request.args.get("value")
    if patient_id not in twin_instances:
        return jsonify({"error": f"Patient {patient_id} not found"}), 404
    
    ## Change the .json file for the patient
    with open(twin_instances[patient_id].param_file, "r") as f:
        params = json.load(f)

        for a,b in params.items():
            for c,d in b.items():
                if c == param:
                    print(f"Changing {param} from {d} to {value}")
                    params[a][c] = value
        
        with open("current.json", "w") as f:
            json.dump(params, f)
    ## Reload the model with the new parameters
    twin_instances[patient_id].model.update_parameters("current.json")
    return jsonify({"status": f"Parameter {param} changed to {value} for patient {patient_id}"})

@app.get('/get_alarms')
def get_alarms():
    """Get the alarms for a specific patient."""
    patient_id = request.args.get("patient_id")
    if patient_id not in twin_instances:
        return jsonify({"error": f"Patient {patient_id} not found"}), 404
    return jsonify(twin_instances[patient_id].model.alarms)

## Call for getting latest data-points for the patient (n=120)
@app.get('/get_latest_data')
def get_latest_data():
    """Get the latest data for a specific patient."""
    patient_id = request.args.get("patient_id")
    if patient_id not in twin_instances:
        return jsonify({"error": f"Patient {patient_id} not found"}), 404
    return jsonify(twin_instances[patient_id].model.data_epoch)

if __name__ == '__main__':
    socketio.run(app, host="0.0.0.0", port=5001)











## Socket.IO event handlers  


@app.route('/patient/<patient_id>')
def view_patient(patient_id):
    """Display monitored values for a specific patient."""
    if patient_id not in twin_instances:
        return f"Patient {patient_id} not found.", 404
    return render_template('patient.html', patient_id=patient_id)

@app.route("/set_resources", methods=["POST"])
def set_resources():
    """Update healthcare resources and setting."""
    resources['residents'] = int(request.form.get("residents", 0))
    resources['nurses'] = int(request.form.get("nurses", 0))
    resources['physicians'] = int(request.form.get("physicians", 0))
    resources['other_staff'] = int(request.form.get("other_staff", 0))
    resources['setting'] = request.form.get("setting", "NICU")  # Default to NICU if not provided

    print("Updated Resources and Setting:", resources)  # Debugging
    return redirect(url_for("home"))

# Socket.IO event handlers
@socketio.on('connect')
def handle_connect():
    print('Client connected')

@socketio.on('disconnect')
def handle_disconnect():
    print('Client disconnected')

@socketio.on('join')
def on_join(data):
    patient_id = data
    join_room(patient_id)
    print(f'Client joined room for patient {patient_id}')


