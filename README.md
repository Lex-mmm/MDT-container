# MDT-container

## Overview

MDT-container is a Flask-based application that simulates digital twin models for patients. The application allows you to start and stop simulations for different patients and view their monitored values in real-time.

## Project Structure




## Installation

1. Clone the repository:
    ```sh
    git clone https://github.com/your-repo/MDT-container.git
    cd ...
    ```

2. Create a virtual environment and activate it:
    ```sh
    python -m venv venv
    source venv/bin/activate  # On Windows use `venv\Scripts\activate`
    ```

3. Install the required packages:
    ```sh
    pip install -r requirements.txt
    ```

## Usage

1. Start the Flask application:
    ```sh
    python main.py
    ```

2. Open your browser and navigate to `http://localhost:5001` to access the home page.

3. Use the form on the home page to start a new simulation for a patient by providing a patient ID and parameter file path.

4. View the running simulations and their monitored values by clicking on the patient links.

## API Endpoints

- **Home Page**
    - `GET /`
    - Displays the list of patients and a form to start a new simulation.

- **Start Simulation**
    - `POST /start_simulation`
    - Starts a simulation for a specific patient.
    - Parameters:
        - `patient_id` (required): The ID of the patient.
        - `param_file` (optional): The path to the parameter file (default: `parameters.json`).

- **Stop Simulation**
    - `POST /stop_simulation/<patient_id>`
    - Stops the simulation for a specific patient.

- **View Patient**
    - `GET /patient/<patient_id>`
    - Displays the monitored values for a specific patient.

## WebSocket Events

- **Connect**
    - Triggered when a client connects to the server.

- **Disconnect**
    - Triggered when a client disconnects from the server.

- **Join**
    - Triggered when a client joins a room specific to a patient.
    - Parameters:
        - `patient_id`: The ID of the patient.

## Digital Twin Model

The `DigitalTwinModel` class in `digital_twin_model_v5.py` is responsible for simulating the patient's digital twin. It uses the parameters from the provided JSON file to initialize the model and runs the simulation in a separate thread.

## License

This project is licensed under the MIT License. See the LICENSE file for details.
