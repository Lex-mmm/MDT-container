let heartRateData = [];
let timeData = [];
const ctx = document.getElementById('heartRateChart').getContext('2d');
const heartRateChart = new Chart(ctx, {
    type: 'line',
    data: {
        labels: timeData,
        datasets: [{
            label: 'Heart Rate',
            data: heartRateData,
            borderColor: 'rgba(75, 192, 192, 1)',
            borderWidth: 1,
            fill: false,
        }]
    },
    options: {
        scales: {
            x: {
                title: {
                    display: true,
                    text: 'Time (s)'
                }
            },
            y: {
                title: {
                    display: true,
                    text: 'Heart Rate (bpm)'
                },
                beginAtZero: true
            }
        }
    }
});

function startDataFetch(patientId) {
    setInterval(() => {
        fetch(`/patient_data/${patientId}`)
            .then(response => response.json())
            .then(data => {
                if (!data.error) {
                    document.getElementById('heart-rate').textContent = data.heart_rate.toFixed(2);
                    document.getElementById('simulation-time').textContent = data.time.toFixed(2);

                    // Update chart data
                    timeData.push(data.time.toFixed(2));
                    heartRateData.push(data.heart_rate.toFixed(2));

                    // Keep the last 50 data points
                    if (timeData.length > 50) {
                        timeData.shift();
                        heartRateData.shift();
                    }

                    heartRateChart.update();
                } else {
                    console.error(data.error);
                }
            })
            .catch(error => console.error('Error fetching data:', error));
    }, 1000); // Fetch data every 1 second
}
