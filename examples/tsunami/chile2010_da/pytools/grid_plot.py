import subprocess


def grid_plot(image_array, outfile):
    from PIL import Image
    new_im = Image.new('RGB',(2400,1900), 'white')
    ens_num = 0
    for i in xrange(0,1800,800):
        for j in xrange(100,1700, 600):
            im = Image.open(image_array[ens_num])
            print im.size
            new_im.paste(im,(i,j))
            ens_num = ens_num + 1

    new_im.save(outfile)

def grid2gif4(image_str, output_gif):
    str1 = 'convert -delay 100 -loop 1 ' + image_str  + ' ' + output_gif
    subprocess.call(str1, shell=True)

def grid2gif5(image_str, output_gif):
    from scitools.std import movie
    movie(image_str,fps=1,output_file=output_gif)
