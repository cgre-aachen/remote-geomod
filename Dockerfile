# install ubuntu with miniconda3
FROM continuumio/miniconda3
# install packages
RUN conda install jupyter numpy matplotlib pandas seaborn theano scikit-image tqdm
# specific version of gdal
RUN conda install -c conda-forge gdal
# steoreonet plots and gempy via pip
RUN pip install mplstereonet

# clone gempy's rgeomod branch
RUN git clone -b rgeomod https://github.com/cgre-aachen/gempy.git

RUN git clone https://github.com/cgre-aachen/remote-geomod.git
