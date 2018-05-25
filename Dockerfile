# install ubuntu with miniconda3
FROM continuumio/miniconda3
# install packages
RUN conda install jupyter numpy matplotlib pandas seaborn theano=1.0.1 scikit-image tqdm nb_conda
# specific version of gdal
RUN conda install -c conda-forge gdal
# steoreonet plots and gempy via pip
RUN pip install mplstereonet
RUN pip install gempy
# clone rgeomod
RUN git clone https://github.com/cgre-aachen/remote-geomod.git
