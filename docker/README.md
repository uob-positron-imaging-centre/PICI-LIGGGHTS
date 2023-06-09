# Docker image of PICI-LIGGGHTS

This Docker image was created to run LIGGGHTS (and the most important Python libraries) in online Jupyter Notebooks - e.g. hosted by Binder.

It is based on the [deal.II](https://hub.docker.com/r/dealii/dealii/) Docker image, which already contains the most complex packages: VTK, MPI, Python, etc. - I have not found a smaller image containing these, so thanks deal.II people!


## Building locally
You can build and test your image locally like this (based on the very helpful [Binder docs](https://mybinder.readthedocs.io/en/latest/tutorials/dockerfile.html))

1. Try building your image.
```bash
docker build -t my-image .
```

2. Try starting a container from the image.
```bash
docker run -it --rm -p 8888:8888 my-image jupyter notebook --NotebookApp.default_url=/lab/ --ip=0.0.0.0 --port=8888
```

3. Inspect the container from terminal.
```bash
docker run -it --rm my-image bash
# what username do i have?
whoami
# what user id do i have?
id -u
# what is the current working directory?
pwd
# who is the owner of the files in the users home directory?
ls -alh ~
```


## Using the Docker image

See [this Binder Environment repository](https://github.com/uob-positron-imaging-centre/ACCES-GranuDrum-Calibration-Env) for an example of using this Docker image to build an environment.

See [this tutorial repository](https://github.com/uob-positron-imaging-centre/ACCES-GranuDrum-Calibration) for interactive examples of using LIGGGHTS in online notebooks.


## For developers

To build the image with a new version of PICI-LIGGGHTS, update the git tag downloaded by `install_liggghts.sh` and the docker tag used in `docker_build_push.sh`, then run `docker_build_push.sh`.