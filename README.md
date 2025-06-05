# MDT-container

## Overview

MDT-container is a Flask-based application designed to simulate digital twin models for patients. This system allows users to initiate, manage, and terminate simulations for individual patients, and to monitor their physiological values in real-time through a web interface.

## Project Structure

```
MDT-container/
├── Alarms/
│   └── AscomMathFunctions.py
├── Visualisation/
│   └── userInterface.py
├── main.py
├── digital_twin_model_v5.py
├── parameters.json
├── requirements.txt
├── Dockerfile
├── DockerStart.sh
└── README.md
```
*(Please adjust the structure above if it's not accurate)*

## Prerequisites

*   Python 3.8+
*   pip (Python package installer)
*   Docker (optional, for containerized deployment)

## Installation

1.  **Clone the repository:**
    ```bash
    git clone <your-repository-url>
    cd MDT-container
    ```

2.  **Create and activate a virtual environment:**
    ```bash
    python3 -m venv venv
    source venv/bin/activate  # On Windows use `venv\Scripts\activate`
    ```

3.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

## Usage

### Running the Application Directly

1.  **Start the Flask application:**
    ```bash
    python main.py
    ```

2.  **Access the application:**
    Open your web browser and navigate to `http://localhost:5001`.

3.  **Using the Interface:**
    *   The home page displays a list of current patients and a form to start new simulations.
    *   To start a new simulation, provide a Patient ID and optionally, a path to a parameter file.
    *   Click on a patient's link to view their running simulation and monitored values.

### Running with Docker

1.  **Build and start the Docker container:**
    Ensure Docker is installed and running on your system.
    ```bash
    ./DockerStart.sh
    ```
    This script should handle building the Docker image and running the container. The application will then be accessible at `http://localhost:5001`.

### Visualizing Results

After running simulations, you can visualize the results using the provided script:

1.  **Run the visualization interface:**
    ```bash
    python Visualisation/userInterface.py
    ```
    Follow the instructions provided by the script to view visualizations.

## API Endpoints

The application provides the following API endpoints:

*   `GET /`: Displays the home page with the list of patients and the simulation form.
*   `POST /start_simulation`: Starts a new simulation.
    *   Parameters: `patient_id` (required), `param_file` (optional).
*   `POST /stop_simulation/<patient_id>`: Stops the simulation for the specified patient.
*   `GET /patient/<patient_id>`: Displays the monitored values for the specified patient.

## WebSocket Events

The application uses WebSockets for real-time updates:

*   **Connect**: Triggered when a client connects.
*   **Disconnect**: Triggered when a client disconnects.
*   **Join**: Triggered when a client joins a patient-specific room.
    *   Parameters: `patient_id`.

## Digital Twin Model

The core simulation logic is handled by the `DigitalTwinModel` class, located in `digital_twin_model_v5.py`. This class takes parameters from a JSON file to initialize and run the patient-specific digital twin simulation in a separate thread.

## Mathematical Functions

The `Alarms/AscomMathFunctions.py` module contains mathematical functions used within the system, particularly for dynamic alarming and data quality checks. These include functions for:
*   Data cleaning (e.g., removing missing values).
*   Quantitative data availability checks.
*   Calculation of mean, median, derivatives, Area Under the Curve (AUC), and integrals.

## License

**Copyright (c) 2025, Dr. LM van Loon**

**All Rights Reserved.**

**Permissions:**
You must obtain explicit prior written permission from Dr. L.M. van Loon (l.m.vanloon@utwente.nl) before any use, copying, modification, merging, publication, distribution, sublicensing, and/or selling of copies of this software and associated documentation files (the "Software").

**Restrictions:**
1.  The Software and its source code **may not be used for training artificial intelligence (AI) models or as input to any AI systems** without explicit prior written permission from Dr. LM van Loon.
2.  Any derivative works or modifications of this Software are subject to the same licensing terms and restrictions, and also require explicit prior written permission.

For any inquiries regarding the use or licensing of this Software, please contact Dr. L.M. van Loon at l.m.vanloon@utwente.nl.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Contact

For questions or support, please contact Dr. LM van Loon at l.m.vanlooN@utwente.nl.