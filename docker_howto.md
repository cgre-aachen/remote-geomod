# Docker for Googlepicks How-to

## Building the Docker image

Building the docker image from the Dockerfile from within the same directory:

    $ docker build . -t <tag>

## Running the Docker image

To get a list of docker images on your system:

    $ docker images

This will provide you with the image-id you need to run the docker image with the following command:

    $ docker run -i -t -p 8888:8888 <image-id>
    
which runs the image with an interactive terminal and port forwarding for Jupyter Notebooks. The image will run inside 
the /app folder, where the whole repository is copied (this can all be modified in the Dockerfile). To now run a 
working, accessible Jupyter Notebook run the following:

    $ jupyter notebook --ip='*' --port=8888 --no-browser --allow-root
    
and use the provided token-link to access the Notebook server inside your browser.

## Handling volumes

    $ docker volume create my-vol
    
