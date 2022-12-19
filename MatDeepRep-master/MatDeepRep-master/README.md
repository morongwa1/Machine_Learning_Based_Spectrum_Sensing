# MatDeepRep: Deep representation learning tool for image classification with Matlab
<img src="https://raw.githubusercontent.com/GKalliatakis/MatDeepRep/master/MatDeepRep%20Overview.png?raw=true" /> </p>


Author: Grigorios Kalliatakis, University of Essex (gkallia@essex.ac.uk)

Copyright 2017, all rights reserved.

Release: v1.0

Licence: BSD (see COPYING file)


### Introduction

----------
For decades, traditional machine learning systems demanded accurate engineering and significant domain expertise in order to design a feature extractor capable of converting raw data (such as the pixel values ofan image) into a convenient internal representation or feature vector from which a classifier could classify or detect patterns in the input. Today, representation learning methods and principally convolutional neural networks (CNNs) [LeCun et al., 1989] are driving advances at a dramatic pace in the computer vision field after enjoying  a  great  success in  large-scale  image recognition  and  object  detection  tasks.	

### Overview

----------
MatDeepRep is a Matlab plugin , built on top of Caffe framework, capable of learning deep representations for image classification using the [BVLC caffe matlab interface (matcaffe)](http://caffe.berkeleyvision.org/tutorial/interfaces.html) & various pretrained .caffemodel binaries.

### Paper and Licencing

----------
Details relating to the supporting paper can be found on my personal research page at:
 * [http://gkalliatakis.com/research/](http://gkalliatakis.com/research/)

If you use this code for your experiments, please cite:

    G. Kalliatakis, S. Ehsan, M. Fasli, A. Leonardis, J. Gall and K. McDonald-Maier
    Detection of Human Rights Violations in Images: Can Convolutional Neural Networks help?
    Computer Vision, Imaging and Computer Graphics Theory and Applications, (VISAPP) Conference, 2017

A copy of the paper is available at:
 * [https://arxiv.org/pdf/1703.04103.pdf](https://arxiv.org/pdf/1703.04103.pdf)

This software is released under the BSD licence (refer to the COPYING file for details). This software is for learning purposes, and not meant for any use in production / commercial purposes. We are always interested in how this software is being used, so if you found this software useful or plan to make a release of code based on or using this package, it would be great to hear from you. Similarly, acknowledgements or citations to the above paper are appreciated.

### Installation

----------
The MatDeepRep function requires the Caffe deep learning framework (available from http://caffe.berkeleyvision.org/), which must be installed and accessible on the MATLAB path (you can follow my step-by-step guide [here](https://github.com/GKalliatakis/Adventures-in-deep-learning/tree/master/Caffe_Installation).

Next, all the latest deep ConvNet models must be downloaded in the required Caffe format from [here](https://docs.google.com/uc?export=download&confirm=WvVf&id=0B98ZKBhlAtp-QUhyaHBnX2NuVU0).
Please note they must be placed inside the default folder of the Caffe installation named 'models'. The safer option will be to extract the downloaded zip file inside the 'models' folder of the Caffe installation. 

In addition, the datasets should be downloaded and stored in a directory of your choice. However, We must pay attention to the structure of the folders which will contain the raw jpg images. We will need 4 different folders: (1) positive training; (2) negative training; (3) positive test and (4) negative test. This happens in order to be able to evaluate both processes of training and testing in the end. The complete structure of the folder with the name 'datasets' (which needs to be inside the matlab/demo folder) is given by the following example:

    datasets/
        MINC2500/
            POS_TRAIN/
                leather/
                    00001.jpg
                    00002.jpg
                    ...

### Usage & examples

----------
Open Matlab and go to Caffe framework. Then add the current folder and subfolders to path. MatDeepRep must be placed and ran from inside the demo folder which can be located inside matlab folder.

The following command will extract fetrures using the ResNet50 model and FMD dataset for the fabric category:

```sh
[code,code_v] = matdeeprep('ReseNet50', 'FMD', 'fabric');
```

You can modify the command according to your needs in order to learn deep features representations for any given dataset.

The final step would be to evaluate your system by feeding the extracted features to a linear SVM classifier, which is beyond the scope of this repository.

### Question and Comments

----------
If you would like to file a bug report or a feature request, use the Github issue tracker.

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)


   [dill]: <https://github.com/joemccann/dillinger>
   [git-repo-url]: <https://github.com/joemccann/dillinger.git>
   [john gruber]: <http://daringfireball.net>
   [@thomasfuchs]: <http://twitter.com/thomasfuchs>
   [df1]: <http://daringfireball.net/projects/markdown/>
   [markdown-it]: <https://github.com/markdown-it/markdown-it>
   [Ace Editor]: <http://ace.ajax.org>
   [node.js]: <http://nodejs.org>
   [Twitter Bootstrap]: <http://twitter.github.com/bootstrap/>
   [keymaster.js]: <https://github.com/madrobby/keymaster>
   [jQuery]: <http://jquery.com>
   [@tjholowaychuk]: <http://twitter.com/tjholowaychuk>
   [express]: <http://expressjs.com>
   [AngularJS]: <http://angularjs.org>
   [Gulp]: <http://gulpjs.com>

   [PlDb]: <https://github.com/joemccann/dillinger/tree/master/plugins/dropbox/README.md>
   [PlGh]:  <https://github.com/joemccann/dillinger/tree/master/plugins/github/README.md>
   [PlGd]: <https://github.com/joemccann/dillinger/tree/master/plugins/googledrive/README.md>
   [PlOd]: <https://github.com/joemccann/dillinger/tree/master/plugins/onedrive/README.md>
