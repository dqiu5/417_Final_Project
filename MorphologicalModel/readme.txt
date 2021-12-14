EECS 417 Project README

HOC Files Needed:
1. init_AssymmetricPulse.hoc (modified for project)
2. makecells.hoc            (modified for project)
3. stim_Asymmetric.hoc      (modified for project)

MOD Files Needed:
1. ipulse1.mod
2. A2hh_k.mod
3. A2hh_na.mod
4. AII_KA.mod
5. AII_KM.mod
6. AII_Na.mod

OTHER Files Needed:
1. swc files in morphology folder 
2. v files in voltage folder
3. 606celltypes.txt         (modified for project)

To Run Simulation:
1. run "mknrndll <modfile name>.mod" to compile the mod files
2. run "./init_AssymmetricPulse.hoc"
3. Open graph->space plot in GUI menu
4. Right click->shape style->show diam
5. Right click->time plot
6. Click on the sections you want to plot
7. In the popped up window, right click->plot what?
8. type in "exIClmp.i" and hit enter
9. Click tools->runcontrol from GUI menu
10. in "Continue till (ms)" type 500 and hit enter
11. Click "Init & Run"

Loizos Paper Citation:
K. Loizos et al, "Increasing electrical stimulation efficacy in degenerated retina: stimulus waveform design in a multiscale computational model," IEEE Transactions on Neural Systems and Rehabilitation Engineering, vol. 26, no. 6, pp. 1111-1120, 2018. Link to the paper



