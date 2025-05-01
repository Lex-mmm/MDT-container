IMAGE_NAME="mdt-pt"
TAG="latest"
DOCKERFILE="dockerfile-mdt-pt"
CONTAINER_NAME="mdt-pt-container"


echo "🚀 Starting container $CONTAINER_NAME..."
docker run -it --name mdt-pt-container mdt-pt:latest