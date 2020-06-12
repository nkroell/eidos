import numpy as np
import imageio
imageio.plugins.freeimage.download()
from cv2 import imwrite, COLOR_BGR2RGB, cvtColor
import os

def read_png_as_zero2one_nparray(path):
    """Reads image from path and returns a scaled numpy array from zero to one.
    
    Args:
        path (Path): Full filename of image.
        
    Returns:
        img as np.ndarray
        
    
    """
    
    img = read_png_as_nparray(path)
    return scale_img_to_01(img)
    

def read_binary_image(path):
    """Reads binary image and returns it as Boolean NumPy array.
    """
    read_image = np.asarray(imageio.imread(path, format="PNG-FI"))
    
    if read_image.dtype == "uint16":
        max_value = 2**16 - 1
    elif read_image.dtype == "uint8":
        max_value = 255
    else:
        max_value = 1
    return read_image[:,:,0] == max_value
    
def read_png_as_nparray(path):
    """Reads 16bit pngs as NumPy arrays.
    
    Args:
        path (Path, String): Path of the image to be loaded.
        
    Returns:
        loaded image as numpy array
    """
    return np.asarray(imageio.imread(path, format="PNG-FI"))


def scale_img_to_01(img):
    """Scales image to 0 to 1 image
    
    Args:
        img (np.ndarray): 16bit or 8bit RGB image
    
    Returns:
        img01 (np.ndarray): image in range 0 to 1.
        
    
    """
    if img.dtype == 'uint16':
        return img/(2**16 - 1)
    elif img.dtype == 'uint8':
        return img/(2**8 - 1)
    else:
        return img

    
def transform_img_2_16bit(img):
    img = img * (2**16 - 1)
    img = img.astype('uint16')
    return img

def transform_img_2_8bit(img):
    img = img * (2**8 - 1)
    img = img.astype('uint8')
    return img

def save_zero_2_one_img_as_16bit(fullfilename_img, img01):
    """ Transforms an numpy array with values from 0 to 1 to uint16 and saves the image in the given path
    
    Args:
        fullfilename_img (Path, String): Full filename of the image (including path).
        img01 (np.ndarray): Image to be saved.
        
    Returns:
        True, if image was saved successfully.
        False, if image saving failed.
    
    """
    img_16bit = transform_img_2_16bit(img01)
    img_16bit = cvtColor(img_16bit, COLOR_BGR2RGB)
    img_saving_successfull = imwrite(str(fullfilename_img), img_16bit)
    return img_saving_successfull

def save_colorpredicited_img(fullfilename_img, img8bit):
    img_saving_successfull = imwrite(str(fullfilename_img), img8bit)
    return img_saving_successfull


def save_zero_2_one_img_as_8bit(fullfilename_img, img01):
    img_8bit = transform_img_2_8bit(img01)
    img_8bit = cvtColor(img_8bit, COLOR_BGR2RGB)
    img_saving_successfull = imwrite(str(fullfilename_img), img_8bit)
    return img_saving_successfull