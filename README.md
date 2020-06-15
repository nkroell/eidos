![alt text](https://github.com/nkroell/eidos/raw/development/images/eidos_logo.png "Eidos Logo")

## Introduction
**Eidos** is the greek word for shape. Quantitative measurement of particle shapes based on binary images is used in many fields, for example chemistry [1], mineral engineering [2], medicine [3] and biology [4]. However, right now there is no comprehensive Python package that allows fast and easy computation of the variety of different shape measurements proposed in the scientific literature (e.g. [5-10]). While there are already powerfull Python packages for image analysis like [scikit-image](https://www.scikit-image.org/) [11] or [opencv](https://www.opencv.org/), combining them means dealing with different coordinate systems, data formats and conventions. Also some shape measurements you find in literature aren't implemented. Here **eidos** comes in.

We unify different coordinate systems, conventions and data formats to calculate - fast and easy - the shape measurements of particles defined as 2D binary images. You insert a 2D binary image into the function and get the shape measurements returned as an Pandas Dataframe.

## How to use eidos?
Using **eidos** is easy - segment your image and insert the binary image `bw` into the function `extract_all_shape_measurements_2d`:

```python

df = extract_all_shape_measurements_2d(bw)

```

As a result you get a [pandas](https://www.pandas.pydata.org/) dataframe where one row representats one found particle and each column an extracted shape measurement.

For details on image segmentation you may look at [scikit-image](https://www.scikit-image.org/). 



## Dependencies

**Eidos** is based on 2D binary images that are defined as NumPy arrays. Where exisisting, we use the powerfull implementation of scikit-image and OpenCV as well as functions given by NumPy and SciPy. **eidos** is developed in Python 3.7.

To use **eidos** you need the following Python packages:
  * [NumPy](https://www.numpy.org/)
  * [SciPy](https://www.scipy.org/)
  * [scikit-image](https://www.scikit-image.org/)
  * [opencv](https://www.opencv.org/)
  * [pandas](https://www.pandas.pydata.org/)

## Literature

[1] Y. M. Lau, N. G. Deen and J. A. M. Kuipers (2013). *Development of an image measurement technique for size distribution in dense bubbly flows*. Chemical Engineering Science, 94, pp. 20–29. https://doi.org/10.1016/j.ces.2013.02.043

[2] T. Andersson, M. J. Thurley and Carlson, J. E. (2012). *A machine vision system for estimation of size distributions by weight of limestone particles*. Minerals Engineering, 25(1), pp. 38–46. https://doi.org/10.1016/j.mineng.2011.10.001

[3] T. M. Nguyen and R. M. Rangayyan (2005). *Shape Analysis of Breast Masses in Mammograms via the Fractal Dimension*. 2005 IEEE Engineering in Medicine and Biology 27th Annual Conference, Shanghai, 2005, pp. 3210-3213. https://doi.org/10.1109/IEMBS.2005.1617159

[4] T. G. Smith, G. D. Lange and W. B. Marks (1996). *Fractal methods and results in cellular morphology — dimensions, lacunarity and multifractals*. Journal of Neuroscience Methods, 69(2), pp. 123–136. https://doi.org/10.1016/s0165-0270(96)00080-5

[5] Deutsches Institut für Normung e. V. (2012). *DIN ISO 9276-6 - Darstellung der Ergebnisse von Par-tikelgrößenanalysen: Teil 6: Deskriptive und quantitative Darstellung der Form und Morphologievon Partikeln*.

[6] M. Pahl, G. Schädel und H. Rumpf (1973). *Zusammenstellung von Teilchenformbeschreibungs-methoden: 1. Teil*. In: Aufbereitungstechnik 14.5, pp. 257–264.

[7] M. Pahl, G. Schädel und H. Rumpf (1973). *Zusammenstellung von Teilchenformbeschreibungs-methoden: 2. Teil*. In: Aufbereitungstechnik 14.10, pp. 672–683.

[8] M. Pahl, G. Schädel und H. Rumpf (1973). *Zusammenstellung von Teilchenformbeschreibungs-methoden: 3. Teil*. In: Aufbereitungstechnik 14.11 , pp. 759–764.

[9] W. Pabst und E. Gregorova (2007). *Characterization of particles and particle systems*.

[10] M. Steuer (2010). *Serial classification*. In: AT Mineral Processing.

[11] S. Walt, J. Schönberger, J. Nunez-Iglesias, F. Boulogne, J. Warner, N. Yager, E. Gouillart, T. Yu and the scikit-image contributors. scikit-image: Image processing in Python. PeerJ 2:e453 (2014) https://doi.org/10.7717/peerj.453

## License
**Eidos** is published under the MIT-License. Your are free to use it. If you have suggestions for improvement or new features feel free to post a new issue.
