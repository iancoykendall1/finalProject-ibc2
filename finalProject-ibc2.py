# -*- coding: utf-8 -*-
"""
math260finalcode.py
Name(s): Ian Coykendall
NetId(s): ibc2
Date: 4/28/2022
"""

#imports
import matplotlib.pyplot as plt

#f functio for u
def fu(u,v,a):
    return (u*(1-u)-a*u*v)

#f function for v
def fv(u,v,B,lam):
    return (lam*v*(1-v)-B*u*v)
    
#RK4 function for population competition model
def RK4(a, B, lam,u0 ,v0, T, N):
    #creates list for t, N+1 long, spaced delta t apart
    tvec = [(i*T)/N for i in range(N+1)]

    #compute delta t
    dt = T/N
    
    #created empty u and v vectors, plugs in initial conditions for 0th entry
    uvec = [0. for i in range(N+1)]
    vvec = [0. for i in range(N+1)]
    uvec[0] = u0
    vvec[0] = v0
    
    #itterates through N and alters u and v vector based on functions and RK4
    for i in range(N):
        #set previous u and v entry conditions      
        un = uvec[i]
        vn = vvec[i]
        
        #computes k1 for u and v
        k1u = dt*fu(un,vn,a)
        k1v = dt*fv(un,vn,B,lam)
        #computes k2 for u and v
        k2u = dt*fu(un+(k1u/2),vn+(k1v/2),a)
        k2v = dt*fv(un+(k1u/2),vn+(k1v/2),B,lam)
        #computes k4 for u and v
        k3u = dt*fu(un+(k2u/2),vn+(k2v/2),a)
        k3v = dt*fv(un+(k2u/2),vn+(k2v/2),B,lam)
        #computes k4 for u and v
        k4u = dt*fu(un+k3u,vn+k3v,a)
        k4v = dt*fv(un+k3u,vn+k3v,B,lam)
        
        #computes new u and v entries based on RK4
        unplus1 = un + (k1u/6) + (k2u/3) + (k3u/3) + (k4u/6)
        vnplus1 = vn + (k1v/6) + (k2v/3) + (k3v/3) + (k4v/6)
        
        #updates i+1 entries of previoulsy empty u and v vectors
        uvec[i+1] = unplus1
        vvec[i+1] = vnplus1
    
    #return u as vector, v as vector, and time values
    return (uvec,vvec,tvec)

"""
main
"""
if __name__ == '__main__':
    #initial conditions given  
    u = [.1,.3,.5,.7,.9]
    v = [.1,.3,.5,.7,.9]
    
    #T and N same as used for hw6
    T = 100
    N = 10000
    
    
    #SCENARIO A
    #inputs given
    aA = .1
    BA = 2.1
    lamA = 2.
    
    #new figure
    plt.figure()
    
    #empty vectors for plotting startpoints
    ustartsA = []
    vstartsA = []
    
    #empty vectors for plotting endpoints
    uendsA = []
    vendsA = []
    
    #iterate through u and v to get each combination of initial conditions
    #plots each on same plot
    for i in range(len(u)):
        for j in range(len(v)):
            #combination of u and v
            u0A = u[i]
            v0A = v[j]
            #plug into RK4
            uvecA,vvecA,tvecA = RK4(aA,BA,lamA,u0A,v0A,T,N)
            #plot u and v output vecotrs, all on same plot
            plt.plot(uvecA,vvecA)
            
            #add start values to start vectors
            ustartsA.append(uvecA[0])
            vstartsA.append(vvecA[0])
            
            #add end values to end vectors
            uendsA.append(uvecA[N])
            vendsA.append(vvecA[N])
    
    #plot end vectors on same plot
    plt.scatter(ustartsA,vstartsA,marker="x",linewidths = 1)
    
    #plot end vectors on same plot
    plt.scatter(uendsA,vendsA,marker="o",edgecolors = 'r',linewidths = 5)
    
    #compute, plot, and print ubar and vbar from formulas
    ubarA = (lamA*(aA-1))/(aA*BA-lamA)
    vbarA = (BA-lamA)/(aA*BA-lamA)
    #not applicable so do not plot
    #plt.scatter(ubarA,vbarA,marker="o",edgecolors = 'b',linewidths = 2)
    print(ubarA,vbarA)
    
    #plot critical points (0,0), (0,1), and (1,0)
    plt.scatter([0,0,1],[0,1,0],marker="o",edgecolors = 'b',linewidths = 2)
    
    #format graph
    plt.title('Scenerio A: Populations v vs u Over Time for Eeach Initial Condition') # Set the title.    plt.xlabel('u') # Set the x-axis label.
    plt.ylabel('v') # Set the y-axis label.
    plt.grid() # add grid
    plt.xticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1]) #set x ticks by .1
    plt.yticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1]) #set y ticks by .1
    plt.savefig('RKscenA.png', bbox_inches='tight') # Save as a png file.
    plt.show() #show all plots for this figure
    
    
    #SCENARIO B
    
    #inputs given
    aB = 1.1
    BB = 1.
    lamB = 2.
    
    #new figure
    plt.figure()
    
    #empty vectors for plotting startpoints
    ustartsB = []
    vstartsB = []
    
    #empty vectors for plotting endpoints
    uendsB = []
    vendsB = []
    
    #iterate through u and v to get each combination of initial conditions
    #plots each on same plot
    for i in range(len(u)):
        for j in range(len(v)):
            #combination of u and v
            u0B = u[i]
            v0B = v[j]
            #plug into RK4
            uvecB,vvecB,tvecB = RK4(aB,BB,lamB,u0B,v0B,T,N)
            #plot u and v output vecotrs, all on same plot
            plt.plot(uvecB,vvecB)
            
            #add start values to start vectors
            ustartsB.append(uvecB[0])
            vstartsB.append(vvecB[0])
            
            #add end values to end vectors
            uendsB.append(uvecB[N])
            vendsB.append(vvecB[N])
    
    #plot end vectors on same plot
    plt.scatter(ustartsB,vstartsB,marker="x",linewidths = 1)
    
    #plot end vectors on same plot
    plt.scatter(uendsB,vendsB,marker="o",edgecolors = 'r',linewidths = 5)
    
    #compute, plot, and print ubar and vbar from formulas
    ubarB = (lamB*(aB-1))/(aB*BB-lamB)
    vbarB = (BB-lamB)/(aB*BB-lamB)
    #not applicable so do not plot
    #plt.scatter(ubarB,vbarB,marker="o",edgecolors = 'b',linewidths = 2)
    print(ubarB,vbarB)
    
    #plot critical points (0,0), (0,1), and (1,0)
    plt.scatter([0,0,1],[0,1,0],marker="o",edgecolors = 'b',linewidths = 2)
    
    #format graph
    plt.title('Scenerio B: Populations v vs u Over Time for Eeach Initial Condition') # Set the title.    plt.xlabel('u') # Set the x-axis label.
    plt.ylabel('v') # Set the y-axis label.
    plt.grid() # add grid
    plt.xticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1]) #set x ticks by .1
    plt.yticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1]) #set y ticks by .1
    plt.savefig('RKscenB.png', bbox_inches='tight') # Save as a png file.
    plt.show() #show all plots for this figure


    #SCENARIO C
    
    #inputs given
    aC = 2.
    BC = 5.
    lamC = 2.
    
    #new figure
    plt.figure()
    
    #empty vectors for plotting startpoints
    ustartsC = []
    vstartsC = []
    
    #empty vectors for plotting endpoints
    uendsC = []
    vendsC = []
    
    #iterate through u and v to get each combination of initial conditions
    #plots each on same plot
    for i in range(len(u)):
        for j in range(len(v)):
            #combination of u and v
            u0C = u[i]
            v0C = v[j]
            #plug into RK4
            uvecC,vvecC,tvecC = RK4(aC,BC,lamC,u0C,v0C,T,N)
            #plot u and v output vecotrs, all on same plot
            plt.plot(uvecC,vvecC)
            
            #add start values to start vectors
            ustartsC.append(uvecC[0])
            vstartsC.append(vvecC[0])
            
            #add end values to end vectors
            uendsC.append(uvecC[N])
            vendsC.append(vvecC[N])
            
    #plot end vectors on same plot
    plt.scatter(ustartsC,vstartsC,marker="x",linewidths = 1)
    
    #plot end vectors on same plot
    plt.scatter(uendsC,vendsC,marker="o",edgecolors = 'r',linewidths = 5)
    
    #add extra plots around critical point (ubar,vbar)
    
    #additional starting points
    u2 = [.1,.2,.3,.4]
    v2 = [.1,.2,.3,.4,.5,.6,.7]
    
    #the rest is the same but for new staritng points...
    #empty vectors for plotting startpoints
    ustartsC2 = []
    vstartsC2 = []
    
    #empty vectors for plotting endpoints
    uendsC2 = []
    vendsC2 = []
    
    #iterate through u and v to get each combination of initial conditions
    #plots each on same plot
    for i in range(len(u2)):
        for j in range(len(v2)):
            #combination of u and v
            u0C2 = u2[i]
            v0C2 = v2[j]
            #plug into RK4
            uvecC2,vvecC2,tvecC2 = RK4(aC,BC,lamC,u0C2,v0C2,T,N)
            #plot u and v output vecotrs, all on same plot
            plt.plot(uvecC2,vvecC2)
            
            #add start values to start vectors
            ustartsC2.append(uvecC2[0])
            vstartsC2.append(vvecC2[0])
            
            #add end values to end vectors
            uendsC2.append(uvecC2[N])
            vendsC2.append(vvecC2[N])
    
    #plot end vectors on same plot
    plt.scatter(ustartsC2,vstartsC2,marker="x",linewidths = 1)
    
    #plot end vectors on same plot
    plt.scatter(uendsC2,vendsC2,marker="o",edgecolors = 'r',linewidths = 5)
    
    #compute, plot, and print ubar and vbar from formulas
    ubarC = (lamC*(aC-1))/(aC*BC-lamC)
    vbarC = (BC-lamC)/(aC*BC-lamC)
    plt.scatter(ubarC,vbarC,marker="o",edgecolors = 'b',linewidths = 2)
    print(ubarC,vbarC)
    
    #plot critical points (0,0), (0,1), and (1,0)
    plt.scatter([0,0,1],[0,1,0],marker="o",edgecolors = 'b',linewidths = 2)
    
    #format graph
    plt.title('Scenerio C: Populations v vs u Over Time for Eeach Initial Condition') # Set the title.    plt.xlabel('u') # Set the x-axis label.
    plt.ylabel('v') # Set the y-axis label.
    plt.grid() # add grid
    plt.xticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1]) #set x ticks by .1
    plt.yticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1]) #set y ticks by .1
    plt.savefig('RKscenC.png', bbox_inches='tight') # Save as a png file.
    plt.show() #show all plots for this figure


    #SCENARIO D
    
    #inputs given
    aD = .5
    BD = 1.
    lamD = 2.
    
    #new figure
    plt.figure()
    
    #empty vectors for plotting startpoints
    ustartsD = []
    vstartsD = []
    
    #empty vectors for plotting endpoints
    uendsD = []
    vendsD = []
    
    #iterate through u and v to get each combination of initial conditions
    #plots each on same plot
    for i in range(len(u)):
        for j in range(len(v)):
            #combination of u and v
            u0D = u[i]
            v0D = v[j]
            u0D2 = u[i]
            v0D2 = v[j]
            
            #plug into RK4
            uvecD,vvecD,tvecD = RK4(aD,BD,lamD,u0D,v0D,T,N)
            #plot u and v output vecotrs, all on same plot
            plt.plot(uvecD,vvecD)
            
            #add start values to start vectors
            ustartsD.append(uvecD[0])
            vstartsD.append(vvecD[0])
            
            #add end values to end vectors
            uendsD.append(uvecD[N])
            vendsD.append(vvecD[N])

    #plot end vectors on same plot
    plt.scatter(ustartsD,vstartsD,marker="x",linewidths = 1)
    
    #plot end vectors on same plot
    plt.scatter(uendsD,vendsD,marker="o",edgecolors = 'r',linewidths = 5)
    
    #compute, plot, and print ubar and vbar from formulas
    ubarD = (lamD*(aD-1))/(aD*BD-lamD)
    vbarD = (BD-lamD)/(aD*BD-lamD)
    plt.scatter(ubarD,vbarD,marker="o",edgecolors = 'b',linewidths = 2)
    print(ubarD,vbarD)
    
    #plot critical points (0,0), (0,1), and (1,0)
    plt.scatter([0,0,1],[0,1,0],marker="o",edgecolors = 'b',linewidths = 2)
    
    #format graph
    plt.title('Scenerio D: Populations v vs u Over Time for Eeach Initial Condition') # Set the title.
    plt.xlabel('u') # Set the x-axis label.
    plt.ylabel('v') # Set the y-axis label.
    plt.grid() # add grid
    plt.xticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1]) #set x ticks by .1
    plt.yticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1]) #set y ticks by .1
    plt.savefig('RKscenD.png', bbox_inches='tight') # Save as a png file.
    plt.show() #show all plots for this figure
    

    #plot u and v vs time for 1 case from each scenerio
    
    #shorter time for better analysis
    T = 30
    
    #A
    #initial condition
    u0 = .1
    v0 = .3
    
    #run RK4
    uvec,vvec,tvec = RK4(aA,BA,lamA,u0,v0,T,N)
    
    #plot
    plt.figure() #new figure
    plt.plot(tvec,uvec,label = 'u') #plot u over time
    plt.plot(tvec,vvec,label = 'v') #plot v over time
    plt.title('Scenerio A Over TIme: u0 = .1, v0 = .3') # Set the title.
    plt.xlabel('u') # Set the x-axis label.
    plt.ylabel('v') # Set the y-axis label.
    plt.legend() #show legend
    plt.savefig('RKscenAtime.png', bbox_inches='tight') # Save as a png file.
    plt.show() #show all plots for this figure
    
    #B
    #initial condition
    u0 = .5
    v0 = .1
    
    #run RK4
    uvec,vvec,tvec = RK4(aB,BB,lamB,u0,v0,T,N)
    
    #plot
    plt.figure() #new figure
    plt.plot(tvec,uvec,label = 'u') #plot u over time
    plt.plot(tvec,vvec,label = 'v') #plot v over time
    plt.title('Scenerio B Over TIme: u0 = .5, v0 = .1') # Set the title.
    plt.xlabel('u') # Set the x-axis label.
    plt.ylabel('v') # Set the y-axis label.
    plt.legend() #show legend
    plt.savefig('RKscenBtime.png', bbox_inches='tight') # Save as a png file.
    plt.show() #show all plots for this figure
    
    #C
    #initial condition
    u0 = .9
    v0 = .7
    
    #run RK4
    uvec,vvec,tvec = RK4(aC,BC,lamC,u0,v0,T,N)
    
    #plot
    plt.figure() #new figure
    plt.plot(tvec,uvec,label = 'u') #plot u over time
    plt.plot(tvec,vvec,label = 'v') #plot v over time
    plt.title('Scenerio C Over TIme: u0 = .9, v0 = .7') # Set the title.
    plt.xlabel('u') # Set the x-axis label.
    plt.ylabel('v') # Set the y-axis label.
    plt.legend() #show legend
    plt.savefig('RKscenCtime.png', bbox_inches='tight') # Save as a png file.
    plt.show() #show all plots for this figure
    
    #D
    #initial condition
    u0 = .7
    v0 = .1
    
    #run RK4
    uvec,vvec,tvec = RK4(aD,BD,lamD,u0,v0,T,N)
    
    #plot
    plt.figure() #new figure
    plt.plot(tvec,uvec,label = 'u') #plot u over time
    plt.plot(tvec,vvec,label = 'v') #plot v over time
    plt.title('Scenerio D Over TIme: u0 = .7, v0 = .1') # Set the title.
    plt.xlabel('u') # Set the x-axis label.
    plt.ylabel('v') # Set the y-axis label.
    plt.legend() #show legend
    plt.savefig('RKscenDtime.png', bbox_inches='tight') # Save as a png file.
    plt.show() #show all plots for this figure