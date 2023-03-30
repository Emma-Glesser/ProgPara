from turtle import width
from PIL import Image
from numpy import block
from numba import cuda
import numba as nb
import numpy as np
import argparse

gauss_matrix = [[1,4,6,4,1],[4,16,24,26,4],[6,24,36,24,6],[4,16,24,26,4],[1,4,6,4,1]]

def CPUversion(filename) :
    # Récupération de l'image et transformation en np_array
    img = Image.open(filename)
    img_array = np.array(img)
    img_array = np.ascontiguousarray(img_array.transpose(1,0,2))
    img_size = img_array.shape
    #Création du np_array du résultat flouté
    img_arraybis = np.zeros(img_size, dtype=np.uint8)
    H = img_array.shape[0]
    L = img_array.shape[1]
    for i in range(H) :
        for j in range(L) :
            R = 0
            G = 0
            B = 0
            for k in range(i-2, i+2) :
                for l in range(j-2, j+2) :
                    # Cas du centre de l'image
                    if (k>=0 and k<H and l>=0 and l<L) :
                        R += img_array[k][l][0] * gauss_matrix[k-(i-2)][l-(j-2)]
                        G += img_array[k][l][1] * gauss_matrix[k-(i-2)][l-(j-2)]
                        B += img_array[k][l][2] * gauss_matrix[k-(i-2)][l-(j-2)]
                    # Cas des bordures
                    else :
                        R += img_array[i][j][0] * gauss_matrix[k-(i-2)][l-(j-2)]
                        G += img_array[i][j][1] * gauss_matrix[k-(i-2)][l-(j-2)]
                        B += img_array[i][j][2] * gauss_matrix[k-(i-2)][l-(j-2)]
            R /= 256
            G /= 256
            B /= 256
            img_arraybis[i][j] = [R,G,B]
    # Transformation de l'array en tableau compatible avec la sauvegarde et sauvegarde
    img_arraybis = img_arraybis.transpose(1,0,2)
    img_arraybis = np.ascontiguousarray(img_arraybis)
    imgbis = Image.fromarray(img_arraybis)
    filename2 = "gaussblur" + filename
    imgbis.save(filename2)

""" 
d_img : image originelle sous forme de tableau 
d_imgbis : image finale à remplir avec les nouvelles valeurs
gauss_kernel : kernel de gauss contenant les coefficient"""
@cuda.jit
def gaussBlurGPU(d_img,d_imgbis,gauss_kernel) :
    H = d_img.shape[0]
    L = d_img.shape[1]
    R = 0
    G = 0
    B = 0
    # Récupération des coordonnées du pixel
    i, j = nb.cuda.grid(2)
    for k in range(i-2, i+2) :
        for l in range(j-2, j+2) :
            if (k>=0 and k<H and l>=0 and l<L) :
                R += d_img[k][l][0] * gauss_kernel[k-(i-2)][l-(j-2)]
                G += d_img[k][l][1] * gauss_kernel[k-(i-2)][l-(j-2)]
                B += d_img[k][l][2] * gauss_kernel[k-(i-2)][l-(j-2)]
            else :
                R += d_img[i][j][0] * gauss_kernel[k-(i-2)][l-(j-2)]
                G += d_img[i][j][1] * gauss_kernel[k-(i-2)][l-(j-2)]
                B += d_img[i][j][2] * gauss_kernel[k-(i-2)][l-(j-2)]
    R /= 256
    G /= 256
    B /= 256
    d_imgbis[i][j][0] = R
    d_imgbis[i][j][1] = G
    d_imgbis[i][j][2] = B

def GPUversion(filename) :
    # Récupération de l'image et tranformation en tableau
    img = Image.open(filename)
    img_array = np.array(img)
    img_array = np.ascontiguousarray(img_array.transpose(1,0,2))
    # Calcul des tailles des threadblocks et gridblocks
    treadblock_size = 1
    grid_size = img.size
    # Création du tableau de l'image finale
    img_bis = np.zeros(img_array.shape,dtype=np.uint8)
    # Copie des tableaux sur le device
    d_img_array = nb.cuda.to_device(img_array)
    d_img_bis = nb.cuda.to_device(img_bis)
    d_gauss = nb.cuda.to_device(gauss_matrix)
    # Appel au kernel
    gaussBlurGPU[grid_size,treadblock_size](d_img_array,d_img_bis,d_gauss)
    cuda.synchronize()
    # Récupération du tableau contenant l'image modifiée
    img_bis = d_img_bis.copy_to_host()
    # Transformation et sauvegarde de la nouvelle image
    img_bis = img_bis.transpose(1,0,2)
    img_bis = np.ascontiguousarray(img_bis)
    imgbis = Image.fromarray(img_bis)
    filename2 = "gaussblurGPU" + filename
    imgbis.save(filename2)

""" 
d_img : image originelle sous forme de tableau 
d_imgbis : image finale à remplir avec les nouvelles valeurs
gauss_kernel : kernel de gauss contenant les coefficient
H et L : hauteur et largeur de l'image
sum et size : somme des coefficient du kernel de gauss et partie entière de la moitié de la taille du kernel"""
@cuda.jit
def gaussBlurGPUGen(d_img,d_imgbis,gauss_kernel,H,L,sum,size) :
    R = 0
    G = 0
    B = 0
    #Récupération des coordonnées du pixel
    i, j = nb.cuda.grid(2)
    for k in range(i-size, i+size) :
        for l in range(j-size, j+size) :
            if (k>=0 and k<H and l>=0 and l<L) :
                R += d_img[k][l][0] * gauss_kernel[k-(i-size)][l-(j-size)]
                G += d_img[k][l][1] * gauss_kernel[k-(i-size)][l-(j-size)]
                B += d_img[k][l][2] * gauss_kernel[k-(i-size)][l-(j-size)]
            else :
                R += d_img[i][j][0] * gauss_kernel[k-(i-size)][l-(j-size)]
                G += d_img[i][j][1] * gauss_kernel[k-(i-size)][l-(j-size)]
                B += d_img[i][j][2] * gauss_kernel[k-(i-size)][l-(j-size)]
    R /= sum
    G /= sum
    B /= sum
    d_imgbis[i][j][0] = R
    d_imgbis[i][j][1] = G
    d_imgbis[i][j][2] = B

def gaussValue(x,y,sigma) :
    exposant = -((x**2 + y**2) / (2*sigma**2))
    return (1 / (2*np.pi*sigma**2)) * np.exp(exposant)

def createGaussKernel(width,sigma) :
    gausskernel = np.zeros((width,width))
    scale = 10000
    center = int(width/2)
    for x in range(width) :
        for y in range(width) :
           gausskernel[x][y] = gaussValue(x-center, y-center, sigma) * scale
    return gausskernel

def sumGaussKernel(gausskernel) :
    sum = 0
    for i in range(gausskernel.shape[0]) :
        for j in range(gausskernel.shape[1]) :
            sum += gausskernel[i][j]
    return sum

def GPUversionGen(inputfile,outputfile,tx,ty,size,sigma) : 
    # Creation du kernel de Gauss en suivant les paramètres
    gausskernel = createGaussKernel(size,sigma)
    sum = sumGaussKernel(gausskernel)
    center = int(size/2)
    # Récupération de l'image et transformation en tableau
    img = Image.open(inputfile)
    img_array = np.array(img)
    img_array = np.ascontiguousarray(img_array.transpose(1,0,2))
    # Assignation des valeurs pour la grille des threads et la grille des blocs
    treadblock_size = (tx,ty,1)
    grid_size = (int(img.size[0]/tx),int(img.size[1]/ty),1)
    # Création de l'image retour
    blurred_img = np.zeros(img_array.shape,dtype=np.uint8)
    # Exportation des tableaux image originale, modifiée et kernel de Gauss sur le device
    d_img_array = nb.cuda.to_device(img_array)
    d_blurred_img = nb.cuda.to_device(blurred_img)
    d_gauss = nb.cuda.to_device(gausskernel)
    # Application de l'opération floutage en parallèle
    gaussBlurGPUGen[grid_size,treadblock_size](d_img_array,d_blurred_img,d_gauss,img_array.shape[0],img_array.shape[1],sum,center)
    cuda.synchronize()
    # Récupération de l'image modifiée et sauvegarde 
    blurred_img = d_blurred_img.copy_to_host()
    blurred_img = blurred_img.transpose(1,0,2)
    blurred_img = np.ascontiguousarray(blurred_img)
    final_img = Image.fromarray(blurred_img)
    final_img.save(outputfile)

def main() :
    nbThreadsx = 16
    nbThreadsy = 16
    gaussKernelSize = 11
    sigma = 100
    inputfile = ""
    outputfile = ""
    parser = argparse.ArgumentParser("Blur a picture")
    parser.add_argument("inputfile", help="the image to blur")
    parser.add_argument("outputfile", help="the modified image")
    parser.add_argument("--tx", type=int, help="number of threads per block on x")
    parser.add_argument("--ty", type=int, help="number of threads per block on y")
    parser.add_argument("--gsize", type=int, help="width of the kernel gauss")
    parser.add_argument("--gsigma", type=int, help="derivation coefficient of the gauss kernel")
    arguments = parser.parse_args()

    inputfile = arguments.inputfile
    outputfile = arguments.outputfile
    if (arguments.tx) :
        nbThreadsx = int(arguments.tx)
    if (arguments.ty) :
        nbThreadsy = int(arguments.ty)
    if (arguments.gsize) :
        gaussKernelSize = int(arguments.gsize)
    if (arguments.gsigma) :
        sigma = int(arguments.gsigma)
    print("Parameters : \ninput :",inputfile,"\noutput :",outputfile,"\nnb threads en x :",nbThreadsx,"\nnb threads en y :",nbThreadsy,"\nsize of kernel :",gaussKernelSize,"\nsigma :",sigma)
    print("Starting modifing the image")
    GPUversionGen(inputfile,outputfile,nbThreadsx,nbThreadsy,gaussKernelSize,sigma)
    print("Done")

main()




