# install ubuntu with miniconda3
FROM continuumio/miniconda3
# install packages
RUN conda install jupyter numpy matplotlib pandas seaborn theano scikit-image
# specific version of gdal
RUN conda install -c conda-forge gdal
# steoreonet plots and gempy via pip
RUN pip install mplstereonet

# make /app the working folder (we can change this)
WORKDIR /app
# copy files from current folder into app to have it all availabe in the image
ADD . /app

# clone gempy's rgeomod branch
RUN git clone -b rgeomod https://github.com/cgre-aachen/gempy.git
