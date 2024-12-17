# main.py

from flask import Flask, request, jsonify, render_template
from flask_socketio import SocketIO, emit, join_room
import threading
from digital_twin_model import DigitalTwinModel  # Import your updated class

app = Flask(__name__, template_folder='templates', static_folder='static')
socketio = SocketIO(app)

# Dictionary to store instances of DigitalTwinModel keyed by patient_id
twin_instances = {}


@app.route('/')
def home():
    """Home page displaying the list of patients."""
    patient_list = list(twin_instances.keys())
    return render_template('index.html', patients=patient_list)


@app.route('/start_simulation', methods=['POST'])
def start_simulation():
    """Start simulation for a specific patient."""
    patient_id = request.form.get("patient_id")
    param_file = request.form.get("param_file", "parameters.json")
    if not patient_id:
        return jsonify({"error": "Missing patient ID"}), 400

    if patient_id not in twin_instances:
        # Create a new DigitalTwinModel instance
        twin_instances[patient_id] = DigitalTwinModel(patient_id, param_file, data_callback=emit_patient_data)

    model = twin_instances[patient_id]
    if not model.running:
        # Start the simulation in a new thread
        threading.Thread(target=model.start_simulation).start()
        return jsonify({"status": f"Simulation started for patient {patient_id}"})
    else:
        return jsonify({"status": f"Simulation already running for patient {patient_id}"})


@app.route('/stop_simulation/<patient_id>', methods=['POST'])
def stop_simulation(patient_id):
    """Stop simulation for a specific patient."""
    if patient_id in twin_instances and twin_instances[patient_id].running:
        twin_instances[patient_id].stop_simulation()
        return jsonify({"status": f"Simulation stopped for patient {patient_id}"})
    else:
        return jsonify({"status": f"No running simulation for patient {patient_id}"})


@app.route('/patient/<patient_id>')
def view_patient(patient_id):
    """Display monitored values for a specific patient."""
    if patient_id not in twin_instances:
        return f"Patient {patient_id} not found.", 404
    return render_template('patient.html', patient_id=patient_id)


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
    print(f"Client joined room: {patient_id}")


# Function to emit data from the DigitalTwinModel to the client
def emit_patient_data(patient_id, data):
    data['patient_id'] = patient_id
    socketio.emit('patient_data', data, room=patient_id)


if __name__ == "__main__":
    socketio.run(app, host="0.0.0.0", port=5001)


# source .venv/bin/activate

