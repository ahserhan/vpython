from vpython import *
# 
# This version is pure Glowscript
#
# -----------------------------------------------
# Function declarations related to visible objects 
# -----------------------------------------------

def set_scene():
    '''
    Set up scene to give up = z-axis, 
    looking down at x-y plane
    '''
    myscene = canvas()
    myscene.height = 500
    myscene.width = 500
    myscene.range = 8 # meters
    myscene.autoscale = False
    myscene.background = color.white
    myscene.forward = vector(-1,-0.5,-0.3)
    myscene.up = vector(0,0,1)
    myscene.align = "left"
    return myscene

def makeplot():
    '''
    Initialize a plot.
    '''
    gd = graph(fast=True, align='left', width=500, height=220, xtitle='Time', ytitle='Amplitude')
    f1 = gcurve(color=color.blue, width=2, visible=show_signal_checkbox.checked)
    # xfield = gcurve(color=color.green, width=1)
    m_perp = gcurve(color=color.red, width=2, visible=False, visible=show_M_perp_checkbox.checked)
    mz = gcurve(color=color.gray(0.5), width=2, visible=show_Mz_checkbox.checked)
    return gd, f1, mz, m_perp

def Clear_plot():
    '''
    Binds to plot clearing button
    '''
    global gd, f1, mz, t, m_perp, time_label
    t = 0.0
    time_label.text='Time: {:.1f}'.format(t)
    f1.delete()
    gd.delete()
    m_perp.delete()
    gd, f1, mz, m_perp = makeplot()

def Stop_plot(b):
    '''
    Binds to plot start/stop button
    '''
    global plotting
    plotting = not plotting
    if plotting: 
        b.text = "Stop plot"
    else: 
        b.text = "Start plot"

def Show_signal(b):
    global f1
    if b.checked:
        f1.visible = True
    else:
        f1.visible = False

def Show_Mz(b):
    global mz
    if b.checked:
        mz.visible = True
    else:
        mz.visible = False

def Show_M_perp(b):
    global m_perp
    if b.checked:
        m_perp.visible = True
    else:
        m_perp.visible = False


def draw_axes(length=20):
    '''
    Draws coodinate axes on scene
    '''
    xaxis = cylinder(pos=vector(-length/2,0,0), axis=vector(length,0,0), 
                     radius=0.05, color=color.gray(0.7) )
    xlabel = label(pos=vector(length/2+0.3,0,0), text='x', box=False )
    yaxis = cylinder(pos=vector(0,-length/2,0), axis=vector(0,length,0), 
                     radius=0.05, color=color.gray(0.7) )
    ylabel = label(pos=vector(0,length/2+0.3,0), text='y', box=False )
    zaxis = cylinder(pos=vector(0,0,-length/2), axis=vector(0,0,length), 
                     radius=0.05, color=color.gray(0.7) )
    zlabel = label(pos=vector(0,0,length/2+0.3), text='z', box=False )

def Run(b):
    '''
    Binds to Pause/Run button
    '''
    global running
    running = not running
    if running: b.text = "Pause"
    else: b.text = "Run"
    
def Reset_M():
    '''
    Binds to reset M button
    '''
    global m
    m = vector(0,0,chi*h.z)

def Reset_time():
    '''
    Binds to reset t button
    '''
    global t, phi, phi1
    Reset_M()
    Clear_plot()
    phi = dphi
    phi1 = dphi1
    
def Reset_all():
    '''
    Resets all values to startup default.  Binds to
    reset all button.
    '''
    set_default()
    

def Pulse_pi_over_2():
    '''
    Binds to pulse pi/2 button
    '''
    global m
    z_temp = (1.0 + sl_h0.value)*h_static_0 - omega()/gamma
    x_temp, y_temp = RF_function( phi1, RFtype )
    m = m.rotate(angle=-angle_pi_over_2, axis=vector(x_temp,y_temp,z_temp))
    # None

def Pulse_pi():
    '''
    Binds to pulse pi button
    '''
    global m
    z_temp = (1.0 + sl_h0.value)*h_static_0 - omega()/gamma
    x_temp, y_temp = RF_function( phi1, RFtype )
    m = m.rotate(angle=-angle_pi, axis=vector(x_temp,y_temp,z_temp))
    # m = m.rotate(angle=-angle_pi, axis=vector(cos(phi),sin(phi),0))
    # None

def Rotate_frame(b):
    '''
    Binds to rotating/lab frame switch button
    '''
    global rotating
    rotating = not rotating
    if rotating: b.text = "<b>Rotating</b> / Lab"
    else: b.text = "Rotating / <b>Lab</b>"

def Pulse_start(b):
    '''
    Binds to RF on/off button
    '''
    global pulsing
    pulsing = not pulsing
    if pulsing: b.text = "<b>RF ON</b>"
    else: b.text = " Cont. "

def RF_select():
    '''
    Binds to Linear/Rotating field menu
    '''
    global RFtype
    RFtype = RF_menu.selected
    
def RF_function(phi, RF_type='rotating'):
    '''
    Queries RF_type and sets value of RF field components
    '''
    if (RF_type == 'rotating'):
        return sl_h1.value*h_static_0*cos(phi), sl_h1.value*h_static_0*sin(phi)
    elif (RF_type == 'linear'):
        return 2.0*sl_h1.value*h_static_0*cos(phi), 0.0
    else:
        return 0.0, 0.0
#
# Set h0 and h1 entry boxes.
def Set_h1(h1_mag):
    '''
    Updates wtext box for h1 field with h1 magnitude
    '''
    wt_h1.text = '{:1.3f}'.format(h1_mag.value)

def Set_h0(h0_mag):
    '''
    Updates wtext box for h0 field with h0 magnitude
    '''
    wt_h0.text = '{:1.3f}'.format(h0_mag.value)

def Set_Delta_omega(domega):
    '''
    Updates wtext box for Delta-omega value
    '''
    wt_dW.text = '{:1.3f}'.format(domega.value)
    
#
# Set T1, T2 in units of 1/f1    
# T_menu = menu( choices=['1.0','4.0','10.0','100.0','1000.0'], bind=T_select )
def T_select():
    '''
    Binds to T_menu for setting T1.  Also updates T2 as % of T1
    '''
    global T,T2
    T = float(T_menu.selected)/log(2.0)
    T2 = T*float(T2_menu.selected)/100.0

def T2_select():
    '''
    Binds to T2_menu to set T2 as fraction of T1.  T2 menu returns % of T1
    '''
    global T2,T
    T2 = T*float(T2_menu.selected)/100.0

        
# -----------------------------------------------
# MAIN PROGRAM EXECUTION
# -----------------------------------------------

scene = set_scene()

run_button = button(text="Run", pos=scene.caption_anchor, bind=Run)
button(text="Reset M", pos=scene.caption_anchor, bind=Reset_M)
# button(text="Reset time", pos=scene.caption_anchor, bind=Reset_time)
button(text="Reset all", pos=scene.caption_anchor, bind=Reset_all)
scene.append_to_caption('\n RF choice: ')
RF_button = button(text=" Cont.  ", pos=scene.caption_anchor, bind=Pulse_start)
button(text="&pi;/2 pulse", pos=scene.caption_anchor, bind=Pulse_pi_over_2)
button(text="&pi; pulse", pos=scene.caption_anchor, bind=Pulse_pi)
scene.append_to_caption('\n Reference frame: ')
frame_button = button(text="Rotating / <b>Lab</b>", pos=scene.caption_anchor, bind=Rotate_frame)

scene.append_to_caption('\n RF <i>H</i><sub>1</sub> amplitude:\n')
sl_h1 = slider(min=0, max=0.4, value=0.05, length=250, bind=Set_h1, left=25, right=15)
wt_h1 = wtext(text='{:1.3f}'.format(sl_h1.value))
scene.append_to_caption(' <i>H</i><sub>0</sub>')

scene.append_to_caption('\n RF detuning &Delta;<i>&omega;</i>:\n')
sl_dW = slider(min=-1.0, max=1.0, value=0.0, length=250, bind=Set_Delta_omega, left=25, right=15)
wt_dW = wtext(text='{:1.3f}'.format(sl_dW.value))
scene.append_to_caption(' &omega;<sub>0</sub>') 

scene.append_to_caption('\n Static field &Delta;<i>H</i><sub>0</sub> shift:\n')
sl_h0 = slider(min=-1.0, max=1.0, value=0.0, length=250, bind=Set_h0, left=25, right=15)
wt_h0 = wtext(text='{:1.3f}'.format(sl_h0.value))
scene.append_to_caption(' <i>H</i><sub>0,static</sub>') 

scene.append_to_caption(' \n <i>H</i><sub>1</sub> type: ')
RF_menu = menu( choices=['rotating','linear'], selected='rotating', bind=RF_select )

scene.append_to_caption(' \n <i>T</i><sub>1</sub> (half-life): ')
T_menu = menu( choices=['4.0','10.0','20.0','40.0','100.0','1000.0','10000.0'], selected='100.0', bind=T_select )

scene.append_to_caption(' \n <i>T</i><sub>2</sub>* (% of <i>T</i><sub>1</sub>): ')
T2_menu = menu( choices=['100','90','50','10','5','1',], selected='100', bind=T2_select )

scene.append_to_caption(' \n\n')

button(text="Clear plot", bind=Clear_plot)
plot_button = button(text="Start plot", bind=Stop_plot)
show_signal_checkbox = checkbox(bind=Show_signal, text='<font color="blue">Pickup coil</font>  ', checked=True)
show_Mz_checkbox = checkbox(bind=Show_Mz, text='<font color="gray">M<sub>z</sub></font>  ', checked=False)
show_M_perp_checkbox = checkbox(bind=Show_M_perp, text='<font color="red">Detector</font>', checked=False)
scene.append_to_caption(' \n\n')

gd, f1, mz, m_perp = makeplot()

# ----------------------------------------------
# CALCULATIONAL FUNCTIONS
# ----------------------------------------------


def omega():
    '''
    Returns RF angular frequency from slider value
    '''
    return omega0*(1.0 + sl_dW.value)

# omega = lambda : omega0*(1.0 + sl_dW.value)

def get_delta_m():
    '''
    Calculate the change in M via the Bloch equation. This is a 2-step process
    Step 1 is to forward shoot to calculate the midpoint of the change.
    Step 2 is to use the estimated value to calculate the new value.
    '''
    # Step 1
    dm_mid = -(gamma*cross(h,m))*0.5*dt  # Cross product part
    m_temp1 = m - chi*h
    dm_mid.z += -(m_temp1.z/T)*0.5*dt # Relaxation part
    dm_mid.y += -(m_temp1.y/T2)*0.5*dt
    dm_mid.x += -(m_temp1.x/T2)*0.5*dt
        
    ## Step 2
    m_temp2 = m+dm_mid
    dm = -(gamma*cross(h,m_temp2))*dt
    m_temp1 = m_temp2 - chi*h
    dm.z += -(m_temp1.z/T)*dt
    dm.y += -(m_temp1.y/T2)*dt
    dm.x += -(m_temp1.x/T2)*dt
    return dm

# -----------------------------------------------
# PARAMETER DECLARATIONS
# -----------------------------------------------
#
# Timestep information
t = 0.0 # start at t=0
ticks_per_s = 64
dt = float(1/ticks_per_s)
#
# These are used to skip visualization frames in
# the plot of RF signal
# fps = 20
# ticks_per_frame = int(ticks_per_s/fps)

# Physical constants
#
chi = 0.8 # static susceptibility
h_static_0 = 5.0 # static magnetic field
h_RF_0 = 0.1*h_static_0 # RF amplitude 
gamma = 2.0 # gyromagnetic ratio
omega0 = gamma*h_static_0 # nominal resonant rotation frequency
omega1 = gamma*h_RF_0 # nominal RF rotation frequency
T0 = float(T_menu.selected)/log(2.0) # relaxation time, in units of 1/(f_1)

# Constant vectors
#
x_hat = vector(1,0,0)
y_hat = vector(0,1,0)
z_hat = vector(0,0,1)
origin_0 = vector(0,0,0)
angle_pi_over_2 = pi/2.0
angle_pi = pi

# 
# Magnetization vector initialization and creation
h0 = vector(0.,0.,h_static_0)
m0 = chi*h0

#
# Loop variable declarations

phi = 0.0 # cumulative phase of m
dphi = 0.0 # increment of phi
phi1 = omega()*t # phase of RF
dphi1 = omega()*dt # increment of phi1
h = h0 # initial field vector
m = m0 # initial magnetization
dm = vector(0.,0.,0.) # initial magnetization increment 


def set_default():
    '''
    Sets default values of all sliders and starting variables
    '''
    global sl_h1, wt_h1, sl_dW, wt_dW, sl_h0, wt_h0
    global RF_menu, T_menu, T2_menu, time_label
    global phi, dphi, phi1, dphi1, h, m, dm
    global running, rotating, pulsing, plotting
    global run_button, RF_button, frame_button
    global scene
    
    running = False
    run_button.text = "Run"
    rotating = False
    frame_button.text = "Rotating / <b>Lab</b>"
    pulsing = False
    RF_button.text = " Cont.  "
    plotting = False
    plot_button.text = "Start plot"
    
    scene.forward = vector(-1,-0.5,-0.3)
    
    sl_h1.value = 0.05
    wt_h1.text ='{:1.3f}'.format(sl_h1.value)
    sl_dW.value = 0.0
    wt_dW.text ='{:1.3f}'.format(sl_dW.value)
    sl_h0.value = 0.0
    wt_h0.text = '{:1.3f}'.format(sl_h0.value)
    RF_menu.selected = 'rotating'
    T_menu.selected = '100.0'
    T2_menu.selected = '100'
    # 
    # Initial values of loop variables
    phi = 0.0 # cumulative phase of m
    dphi = 0.0 # increment of phi
    phi1 = omega()*t # phase of RF
    dphi1 = omega()*dt # increment of phi1
    h = h0 # initial field vector
    m = m0
    dm = vector(0.,0.,0.) # initial magnetization increment 
    T_select()
    Clear_plot()
    time_label.text='Time: {:.1f}'.format(t)


# 
# Initial control states
running = False
rotating = False
pulsing = False
plotting=False
T = T0 # relaxation time (T1 or T2)
T2 = T
DH0 = 0.0 # static field sweep

# --------------------------------------------------
# MAIN PROGRAM EXECUTION
# --------------------------------------------------

# Draw coordinate axes
#
draw_axes(length=11)

label(pos=vec(10,35,0), align='left', pixel_pos=True, 
      text='Green: Magnetic field <b><i>B</i></b>', 
      color=vector(0,0.6,0), opacity=0, box=False, height=18)
      
label(pos=vec(10,15,0), align='left', pixel_pos=True, 
      text='Blue: Net magnetization <b><i>M</i></b>', 
      color=vector(0,0,0.9), opacity=0, box=False, height=18)
      
time_label = label(text='Time: {:.1f}'.format(t), pixel_pos=True, box=False,
                   align='left', opacity=0, height=18, pos=vector(300,20,0))

# Draw magnetization and field vectors
#
M = arrow(pos=origin_0, axis=m0, shaftwidth=0.3, round=True, color=color.blue)
B = arrow(pos=origin_0, axis=h0, shaftwidth=0.3, round=True, color=color.green)


# Main loop
while True:
    rate(ticks_per_s)
    
    ## Set value of h in the current step
    h.z = (1.0 + sl_h0.value)*h_static_0
    
    ## Add RF component, if needed
    if pulsing:
        h.x, h.y = RF_function( phi1, RFtype )
    else:
        h.x = h.y = 0.0
    
    ## Set the axis (length and direction) of the B vector 
    ## Adjust the B (field) vector axis to effective value when in the rotating frame
    if rotating:
        B.axis = h - vector(0,0,omega()/gamma)
    else:
        B.axis = h
    
    ## Set the axis of the M (magnetization) vector to the current value
    M.axis = m
   
    if running: # Update the vectors by recalculating their axes
        if rotating: # Update the camera position 
            scene.forward = scene.forward.rotate(angle=dphi1, axis=z_hat)

        ## Get the change in m
        # dm = get_delta_m()
        
        ## Update all steps
        # m += dm
        dphi = -gamma*h.mag*dt
        m = m.rotate(angle=dphi, axis=h)
        m.x -= ((m - chi*h).x/T2)*dt
        m.y -= ((m - chi*h).y/T2)*dt
        m.z -= ((m - chi*h).z/T)*dt
        phi += dphi
        dphi1 = -omega()*dt
        phi1 += dphi1
        t += dt
        time_label.text='Time: {:.1f}'.format(t)
        m_signal = m.y
        m_transverse = sqrt(m.x**2 + m.y**2)
        
        ## If plotting, update the curve
        if plotting: 
            f1.plot([t, m_signal])
            mz.plot([t, m.z])
            m_perp.plot([t, m_transverse])
