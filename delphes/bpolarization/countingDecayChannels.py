import ROOT, math
from array import array
from collections import defaultdict
execfile("../loadDelphes.py")

def getDecayChannels(*args):
    if(len(args)<3):
        print 'Please give me filedir and filename and pid'
        return
    filedir=args[0]
    filename=args[1]
    pid=args[2]
    
    tfiles = ROOT.TFile(filedir+filename+".root")
    trees=tfiles.Get("Delphes")
    
    #without sign
    channels=defaultdict()
    eventnumber=0
    particlenumber=0
    decaynumber=0
    for iev, event in enumerate(trees):
        eventnumber+=1
        for p in event.Particle:
            if p.Status>30: continue
            if abs(p.PID)!=pid: continue
            sign=p.PID/abs(p.PID)
            particlenumber+=1
            daus=[]
            if p.D1!=-1:
                if p.D2!=-1:
                    for i in range(p.D1,p.D2+1):
                        if event.Particle[i].Status < 30:
                            daus.append(event.Particle[i].PID)
                elif p.D2==-1: daus.append(event.Particle[p.D1].PID)
            daus.sort(key=lambda x:abs(x))
            daus=tuple(daus)
            if channels.has_key(daus):
                channels[daus]+=1
            elif not channels.has_key(daus):
                channels[daus]=1.0
    report=''
    report+='From '+str(eventnumber)+' events \n'
    report+='Particle '+str(pid)+' decay report \n'
    report+='# of Particle '+str(pid)+' : ' +str(particlenumber)+'\n'
    
    keys=channels.keys()
    keys.sort(key=lambda x:abs(x[0]))
    keys.sort(key=lambda x:len(x))
    
    for k in keys:
        report+=str(k)+" : "+str(int(channels[k]))+', fraction('+str(int(channels[k]/particlenumber*100))+')\n'
    print report
    f=open('DecayChannelReport.txt','w')
    f.write(report)
    f.close()

    return channels

filedir = "/mnt/delphes/"
filename = 'tsWbW_100k'
#filename = 'tsWbW_20k'
#filename = 'tsWbW_1k'
pid = 5122
channels=getDecayChannels(filedir,filename,pid)
