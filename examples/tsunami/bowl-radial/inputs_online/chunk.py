import numpy as np

def chunk_write(forlist,type1="init"):
    from shutil import copyfile
    #for i in range(1,num_ens+1):
    for i in forlist:
        if type1=="init":
            srcfile = "ens_"+str(i)+".txt"
            destfile = "ens_"+str(i)+"_mod.txt"
            print "Saving num_ens = ", i

        if type1=="obs":
            srcfile = "obs_step"+str(i)+".txt"
            destfile = "obs_step"+str(i)+"_mod.txt"
            print "Saving obs_step = ", i

        if type1=="testing":
            srcfile = "masked_topo.txt"
            destfile = "total_array2"

        copyfile(srcfile,destfile)
        z = np.loadtxt(srcfile)

        z1 = np.vsplit(z,2)
        z2 = np.split(z1[0],2,axis=1)
        z3 = np.split(z1[1],2,axis=1)
        chunk1 = np.ravel(z2[0])
        chunk2 = np.ravel(z2[1])
        chunk3 = np.ravel(z3[0])
        chunk4 = np.ravel(z3[1])
        stacked_array = np.row_stack((chunk1.T, chunk2.T, chunk3.T, chunk4.T))
        stacked_array = np.ravel(stacked_array)
        np.savetxt(srcfile,stacked_array)
#    with open(srcfile,"ab") as f:
#        np.savetxt(f,chunk1, fmt='%12.10f')
#        np.savetxt(f,chunk2,fmt='%12.10f')
#        np.savetxt(f,chunk3,fmt='%12.10f')
#        np.savetxt(f,chunk4,fmt='%12.10f')


if __name__=="__main__":
    chunk_write([1],type1="testing")
