# install ubuntu with miniconda3
FROM continuumio/miniconda3
# install packages
RUN conda install jupyter numpy matplotlib pandas seaborn theano scikit-image
# install special (working) build of vtk 7 for py3
RUN conda install -c clinicalgraphics vtk
RUN conda install -c conda-forge gdal
# steoreonet plots and gempy via pip
RUN pip install pymc mplstereonet gempy

# make /app the working folder (we can change this)
WORKDIR /app
# copy files from current folder into app to have it all availabe in the image
ADD . /app

# command to run
# docker run -i -t -p 8888:8888 <IMAGE>
# jupyter notebook --ip='*' --port=8888 --no-browser --allow-root
