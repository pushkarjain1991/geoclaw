from PIL import Image
import subprocess


for folder_num in range(9):
    filename = '_output/_output_' + str(folder_num) + '_for/'
    shell_str = 'python /workspace/Pushkar/clawpack/visclaw/src/python/visclaw/plotclaw.py _output/_output_' + str(folder_num) + '_for _plots_' + str(folder_num) + ' ./setplot.py'
    subprocess.call(shell_str, shell=True)

#for t in range(1,26,1):
for t in range(3,24,1):
    new_im = Image.new('RGB',(3000,3000))
    counter = 0
    for i in xrange(0,3000,1000):
        for j in xrange(0,3000,1000):
            im = Image.open('_plots_' + str(counter) + '/frame' + str(t).zfill(4) +'fig0.png')
            #im.thumbnail(300,300)
            new_im.paste(im,(i,j))
            counter =  counter + 1

    new_im.save('grid_img_' + str(t) + '.png')
