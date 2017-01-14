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

def grid2gif2(image_array, output_gif):
    import imageio
    
    with imageio.get_writer(output_gif, mode='I', duration=0.5) as writer:
        for filename in image_array:
            image = imageio.imread(filename)
            writer.append_data(image)
