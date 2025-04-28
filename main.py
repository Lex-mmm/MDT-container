from flask import Flask, request, jsonify, render_template, redirect, url_for
from flask_socketio import SocketIO, emit, join_room
import threading
from digital_twin_model import DigitalTwinModel  # Import your updated class

from Inference.inference_calc import Inference  # Import inference-generation class
import json

app = Flask(__name__, template_folder='templates', static_folder='static')

# Dictionary to store instances of DigitalTwinModel keyed by patient_id
twin_instances = {}

## Define global class for patient-level information
class Patient:

    def __init__(self, patient_id, param_file="parameters.json", pat_char=None, sleep=True):

        ## Initialize starting instance variables
        self.patient_id = patient_id
        self.param_file = param_file
        self.running = False
        self.data_epoch = []
        self.model = DigitalTwinModel(patient_id, param_file, data_callback=None, 
                                     sleep=sleep)

    def startup(self):
        """Start the simulation."""
        self.running = True
        self.model.start_simulation()


    def stop_patient(self):
        """Stop the simulation."""
        self.running = False
        self.model.stop_simulation()

    def add_disease(self, disease, disease_severity):
        """Add a disease to the patient."""
        self.model.add_disease(disease, disease_severity)









#### API ROUTING ENDPOINTS ####


@app.route('/start_simulation', methods=['POST'])
def start_simulation():
    """Start simulation for a specific patient."""
    ## Allow for back-end submission of patient_id and param_file
    patient_id = request.args.get("patient_id")
    param_file = request.args.get("param_file")
    sleep = request.args.get("sleep")

    if not patient_id:
        return jsonify({"error": "Missing patient ID"}), 400

    if patient_id not in twin_instances:
        # Create a new patient instance
        twin_instances[patient_id] = Patient(patient_id= patient_id, param_file= param_file, sleep=sleep)
    

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
    ## Update the parameter value in the model
    ## PM: Checks for viable parameter ranges
    
    output = twin_instances[patient_id].model.update_param(patient_id, param, value)
    if 'Error' in output:
        return jsonify({"error": output})
    else: 
        return jsonify({"status": f"Parameter {param} changed to {value} for patient {patient_id}"})

## Call for getting alarms for the patient (last n=60)
@app.get('/get_alarms')
def get_alarms():
    """Get the alarms for a specific patient."""
    patient_id = request.args.get("patient_id")
    if patient_id not in twin_instances:
        return jsonify({"error": f"Patient {patient_id} not found"}), 404
    return jsonify(twin_instances[patient_id].model.alarms)

@app.route('/infer_swap', methods = ['POST'])
def infer_swap():
    patient_id = request.args.get("patient_id")
    if patient_id not in twin_instances:
        return jsonify({"error": f"Patient {patient_id} not found"}), 404
    inference_level = request.args.get("inference_level")
    param = request.args.get("param")
    twin_instances[patient_id].model.inference.change_inference(param, inference_level)
    return jsonify({"status": f"Inference level for {param} changed to {inference_level} for patient {patient_id}"})


## Call for getting latest data-points for the patient (n=120)
@app.get('/get_latest_data')
def get_latest_data():
    """Get the latest data for a specific patient."""
    patient_id = request.args.get("patient_id")
    if patient_id not in twin_instances:
        return jsonify({"error": f"Patient {patient_id} not found"}), 404
    return jsonify(twin_instances[patient_id].model.data_epoch)

if __name__ == '__main__':
    print("Starting server")
    app.run(debug=True, host='127.0.0.1', port=5001)

    

