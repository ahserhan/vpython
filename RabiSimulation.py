from vpython import *

#
# -----------------------------------------------
# PARAMETER DECLARATIONS
# -----------------------------------------------
#

# Timestep information
ticks_per_s = 64
pct_speed = 1.0
default_dt = float(1/ticks_per_s)
dt = default_dt 
current_dt = default_dt
t = 0.0 # Clock variable
t_plot = 0.0  # used for plotting and timing in "instant" pulse mode

# Physical constants
chi = 0.8 # static susceptibility
omega0 = 4*pi # f_0 = 2 Hz
gamma = 2.0 # Note: H_0 = omega0/gamma is static magnetic field

# initialization defaults
FWHM = sqrt(2.0*log(2.0))/1.0 # spread of static field, as fraction of field
delta = 1.20 # stretch of pi pulse for M-G effect
T0 = 1.0e6 # relaxation time, in terms of relaxation half-life, here "infinite".
Nspins = 3000 # number of spins
 
# Field vector initialization and creation
Om0 = vector(0.,0.,omega0)
m0 = chi*Om0/gamma

# List objects to include spin vectors, local field values for each spin.
spins = []   # the individual spin vectors
omega_local = [] # the local field for each spin
lil_ems = [] # the sphere/arrow objects for each spin
Tpulse = [] # list of pulse start times

# -----------------------------------------------
# FUNCTION DEFINITIONS
# -----------------------------------------------

def set_scene(size=530):
    '''
    Set up scene to give up = z-axis, 
    looking down at x-y plane
    '''
    myscene = canvas()
    myscene.align = 'left'
    myscene.height = size # pixels
    myscene.width = size
    myscene.range = 8.5 # units of field
    myscene.autoscale = False
    myscene.background = color.white
    myscene.forward = vector(-1,-0.5,-0.3)
    myscene.up = vector(0,0,1)
    return myscene

def draw_axes(length=13):
    '''
    Draws coodinate axes on scene and bundles into compound object
    '''
    xaxis = cylinder(pos=vector(-length/2,0,0), axis=vector(length,0,0), 
                     radius=0.05, color=color.gray(0.7) )
    xlabel = text(pos=-xaxis.pos+vector(0,0,0.2), text="x'", 
                  axis=vector(-1,0,0), align='center', up=vector(0,0,1),
                  height=0.5, color=color.gray(0.8))
    yaxis = cylinder(pos=vector(0,-length/2,0), axis=vector(0,length,0), 
                     radius=0.05, color=color.gray(0.7) )
    ylabel = text(pos=-yaxis.pos+vector(0,0,0.2), text="y'", 
                  axis=vector(0,1,0), align='center', up=vector(0,0,1),
                  height=0.5, color=color.gray(0.8))
    zaxis = cylinder(pos=vector(0,0,-length/2), axis=vector(0,0,length), 
                     radius=0.05, color=color.gray(0.7) )
    zlabel = text(pos=-zaxis.pos+vector(0,0,0.2), align='center', text='z', 
                  axis=vector(0,1,0), up=vector(0,0,1), height=0.5, 
                  color=color.gray(0.8))
    return compound([xaxis, yaxis, zaxis, xlabel, ylabel, zlabel],
                    axis=vector(1,0,0), origin=vector(0,0,0))

def makeplot(size=510):
    '''
    Initialize the plot area and arrays to plot
    '''
    gd = graph(fast=True, align='left', width=size, height=0.5*size,
               xtitle='Time', ytitle='Amplitude')
    f1 = gcurve(color=color.blue, width=2, 
                visible=show_signal_checkbox.checked)
    Mx = gcurve(color=color.red, width=2,
                    visible=show_Mx_checkbox.checked)
    My = gcurve(color=color.green, width=2,
                    visible=show_My_checkbox.checked)
    mz = gcurve(color=color.gray(0.5), width=2,
                visible=show_Mz_checkbox.checked)
    mr_ave = gcurve(color=color.orange, width=2,
                    visible=show_Mr_ave_checkbox.checked)
    return gd, f1, mz, Mx, My, mr_ave
    

def spin_initialize(Hwidth=0.04, N=300): 
    '''
    Initialize the set of spins with their objects, object color
    and local value of the z-field
    '''
    global spins, omega_local, lil_ems
    
    # Clear out the global variables
    for em in lil_ems:
        em.visible = False
        del em
    spins.clear()
    omega_local.clear()
    lil_ems.clear()

    # Set the range array that contains the variation in the z-field
    field_range =[]  
    field_range.append(0.0)
    for i in range(N-1):
        field_range.append(Hwidth*normal_random())

    # This sets up the color variation so that the "fastest" spins
    # are colored red, and the "slowest" are dark blue.    
    max_field_low = abs(min(field_range))
    max_field_hi = abs(max(field_range))
    max_field = max(max_field_hi, max_field_low)
    if max_field < 1e-6: max_field = 1.0  # Don't blow up if max field is zero.

    # Fill the global lists spins, omega_local and lil_ems
    for i in range(N):
        local_field = vector(0.0, 0.0, field_range[i])
        omega_local.append(local_field)
        local_color = 0.55*(1.0 - abs(field_range[i])/max_field)
        if is_arrow:
            lil_ems.append(arrow(pos=vector(0,0,0), axis=m0, shaftwidth=0.02, 
            round=True, color=color.hsv_to_rgb(vector(local_color,1,1))))
        else:
            lil_ems.append(simple_sphere(pos=m0, radius=0.12, 
            color=color.hsv_to_rgb(vector(local_color,1,1))))
        spins.append(m0)                 
    return N


# *** Mainly mathematical functions ***

def Optimal(hx=omega0, hy=0.0, hz=0.0):
    '''
    "Optimal" pulse length gives the largest FID regardles of the detuning
    frequency.  Equivalent to a pi/2 pulse on resonance.
    '''
    cot_theta_sq = hz**2/(hx**2 + hy**2)
    Omega_mag = sqrt(hx**2 + hy**2 + hz**2)
    if (cot_theta_sq > 1):
        return pi/Omega_mag
    else:
        return (pi - acos(cot_theta_sq))/Omega_mag

def Pulse_length(Pulse_type='pi'):
    '''
    Calculate pulse length.  Accurate only for steady fields,
    approximate for oscillating field.
    '''
    global delta
    
    Omx, Omy = RF_function(0,'x-axis')
    Omz = Domega()
    Omega = vector(Omx,Omy,Omz)
    if (Pulse_type == 'pi/2'):
        return pi/(2.0*Omega.mag)
    elif (Pulse_type == 'pi'):
        return pi/(Omega.mag)
    elif (Pulse_type == 'pi+error'):
        return (pi*delta)/(Omega.mag)
    else:
        return Optimal(Omx, Omy, Omz)    

def RF_function(phi, RF_type='x-axis'):
    '''
    RF function along x' or y' axes, or linearly oscillating 
    along x (lab frame) axis.  Binds to RF amplitude slider.
    '''
    if (RF_type == 'x-axis'):
        return sl_h1.value*omega0, 0.0
    elif (RF_type == 'oscillating'):
        return sl_h1.value*omega0*(1+cos(2.0*phi)), sl_h1.value*omega0*sin(2.0*phi) 
    elif (RF_type == 'y-axis'):
        return 0.0, sl_h1.value*omega0
    else:
        return 0.0, 0.0

def ave_vector(v_list=[m0]):
    '''
    Calculate average magnetization vector
    '''
    ave_v = vector(0.,0.,0.)
    ave_v_perp = 0.0
    for v in v_list:
        ave_v = ave_v + v
        ave_v_perp = ave_v_perp + sqrt(v.x**2+v.y**2)
    return ave_v/len(v_list), ave_v_perp/len(v_list)

def Domega():
    '''
    This makes the angular frequency of the RF proportional to both 
    H1 magnititude and the detuning. It allows the pulse rotation vector
    direction to be independent of the rate of rotation.
    '''
    return sl_h1.value*omega0*-sl_dW.value

def normal_random():
    '''
    Box-Muller transform to get normal distribution from uniform 
    0-1 random distribution
    '''
    u1 = random()
    u2 = random()
    return sqrt(-2.0*log(u1))*cos(2.0*pi*u2)
    
# *** Mainly functions bound to controls ***

# Row 1 - Initialize spins ------------------------

def Set_spin_shape(b):
    '''
    Select object type to represet spins.
    Binds to Balls/Arrows button.
    '''
    global prog_ctl, is_arrow
    
    prog_ctl = 'init'
    is_arrow = not is_arrow
    if is_arrow: b.text = "Balls/<b>Arrows</b>"
    else: b.text = "<b>Balls</b>/Arrows"
        
def Spin_select():
    '''
    Sets the number of spins to create
    Binds to Spin_number menu.
    '''
    global prog_ctl, newspins
    
    prog_ctl = 'init'
    newspins = int(Spin_number.selected)
    
def FWHM_select():
    '''
    Sets T2* width.
    Binds to T2* entry box fwhm_w
    '''
    global prog_ctl, FWHM
    
    prog_ctl = 'init'
    FWHM = (sqrt(2.0*log(2.0)))/float(fwhm_w.text)

# Row 2/3 - Program Control ------------------------- 

def Run(b):
    '''
    Binds to the Pause/Run button
    '''
    global running
    
    running = not running
    if running: b.text = "Pause"
    else: b.text = "Run"
    
def Reset_M():
    '''
    Binds to Reset M button.  Also clears any pending pulses
    in a pulse sequence.
    '''
    global spins, Npulse, pulse_time, wt_next_pulse
    
    Npulse = 0
    pulse_time = 0.0
    wt_next_pulse.text = ' ' 
    for i in range(len(spins)):
        spins[i] = m0

def Clear_plot():
    '''
    Binds to Reset t button.  Resets time to zero and
    reinitializes all of the plotting variables
    '''
    global gd, f1, Mx, My, t_plot, t, mz, mr_ave 
    global scene, pulse_time, wt_time, tip
    
    # scene.forward = vector(-1,-0.5,-0.3)
    # print_options(clear=True)
    if spin_0_trail: tip.clear_trail()
    t = 0.0 #dt
    pulse_time = 0.0
    t_plot = 0.0 # dt
    wt_time.text ='Time: {:.1f}'.format(t_plot)
    f1.delete()
    mz.delete()
    Mx.delete()
    My.delete()
    mr_ave.delete()
    gd.delete()
    gd, f1, mz, xfield, yfield, mr_ave = makeplot()
    
def Rotate_frame(b):
    '''
    Binds to rotating/lab frame switch button
    '''
    global in_rotating_frame, RF_frame_checkbox
    
    in_rotating_frame = not in_rotating_frame
    if in_rotating_frame: 
        b.text = "<b>Rotating</b>/Lab"
        RF_frame_checkbox.disabled = False
    else: 
        b.text = "Rotating/<b>Lab</b>"
        RF_frame_checkbox.disabled = True

def Set_sim_speed(sl_sim_speed):
    '''
    Binds to Simulation speed slider
    '''
    global pct_speed, current_dt
    
    pct_speed = 0.01*int(sl_sim_speed.value)
    # current_dt = default_dt*pct_speed
    wt_sim_speed.text = '{:d}%'.format(int(sl_sim_speed.value))

def Rotating_frame_select(b):
    '''
    When checked, rotating frame is bound to RF frequency only.
    When unchecked, rotating frame switches between free-field
    and RF frequency, depending on whether RF is on.
    '''
    global RF_frame_only
    
    if b.checked:
        RF_frame_only = True
    else:   
        RF_frame_only = False
        

# Row 4 - Continuous RF

def RF_start(b):
    '''
    Binds to RF on/off button
    '''
    global RF_on
    
    RF_on = not RF_on
    if RF_on: b.text = "<b>ON</b>/Off"
    else: b.text = "On/<b>OFF</b>"

def RF_select():
    '''
    Binds to RF type menu. Chooses which
    type of RF: along x', along y' or oscillating
    along lab frame  (nonrotating) x axis.
    Used with RF_function().
    '''
    global RFtype
    
    RFtype = RF_menu.selected

# Set h1 and delta-omega entry boxes.
def Set_h1(h1_mag):
    wt_h1.text = '{:1.3f}'.format(h1_mag.value)
    wt_dW_0.text = '{:1.3f}'.format(sl_dW.value*h1_mag.value)

def Set_Delta_omega(domega):
    wt_dW.text = '{:1.3f}'.format(domega.value)
    wt_dW_0.text = '{:1.3f}'.format(domega.value*sl_h1.value)

def T_select():
    '''
    Binds to T_menu for setting T1.  Also updates T2 as % of T1
    '''
    global T,T2
    if (T_menu.selected == 'Infinite'):
        T = 1.0e6/log(2.0)
    else: 
        T = float(T_menu.selected)/log(2.0)
    T2 = T*float(T2_menu.selected)*0.010 

def T2_select():
    '''
    Binds to T2_menu to set T2 as fraction of T1.  T2 menu returns % of T1
    '''
    global T2
    T2 = T*float(T2_menu.selected)*0.010


def Interval_select():
    global sequence_interval
    sequence_interval = float(Interval_menu.selected)
    # print('Pulse interval:',sequence_interval)
    
def Apulse_select():
    global Apulse_type
    Apulse_type =(Apulse_menu.selected, A_RF_menu.selected)
    # print('A Pulse:', Apulse_type)

def Bpulse_select():
    global Bpulse_type
    Bpulse_type = (Bpulse_menu.selected,B_RF_menu.selected)
    # print('B Pulse:', Bpulse_type)

# The number of pulses being greater than 0 actually starts the 
# pulse sequence, so this is just a dummy function    
def NumB_pulse_select():
    pass
   
def Sequence_start():
    global Npulse, Tpulse, sequence_state
    sequence_state = 'Fire_A'
    Npulse = 1 + int(NumB_pulse_menu.selected)
    Tpulse = [t_plot]
    for i in range(Npulse-1):
        if i == 0:
            nextpulse = sequence_interval + Tpulse[i]
        else:
            nextpulse = 2*sequence_interval + Tpulse[i]
        Tpulse.append(nextpulse)
    # print('Npulse: ', Npulse)
    # print('pulse times: ', Tpulse)

# "Instant" means that during the pulse, the spins are prevented from evolving
# according to their local field.  An Instant pulse gives the same result
# regardless of the size of RF amplitude.
def Instant_select():
    global instant_pulse
    instant_pulse = Instant_checkbox.checked
    
def Apulse_on():
    global Apulse_yes
    Apulse_yes = Apulse_checkbox.checked
    
def Bpulse_on():
    global Bpulse_yes
    Bpulse_yes = Bpulse_checkbox.checked


# Stop/Start plotting    
def Stop_plot(b):
    global plotting
    plotting = not plotting
    if plotting: 
        b.text = "Stop plot"
    else: b.text = "Start plot"

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

def Show_Mx(b):
    global Mx
    if b.checked:
        Mx.visible = True
    else:
        Mx.visible = False

def Show_My(b):
    global My
    if b.checked:
        My.visible = True
    else:
        My.visible = False
    
def Show_Mr_ave(b):
    global mr_ave
    mr_ave.visible = b.checked

def delta_select():
    '''
    Binds to w_delta to set stretch in pi pulse
    '''
    global delta
    delta = 1.0 + 0.01*float(w_delta.text)


# -----------------------------------------------
# MAIN PROGRAM EXECUTION
# -----------------------------------------------
# Initialize the plot

# Set scene first.  Necessary in order to put controls to right of scene
scene = set_scene()

# Row 1 -----------------------------------------
scene.append_to_caption('<hr width="530" align="left">')
scene.append_to_caption('  INITIALIZE SPINS - Shape:')
button(text="<b>Balls</b>/Arrows", pos=scene.caption_anchor, bind=Set_spin_shape)
scene.append_to_caption(' No.:')
Spin_number = menu( choices=['1','5','10','100','300','1000','3000'], selected='3000', 
                    pos=scene.caption_anchor, bind=Spin_select )
scene.append_to_caption(' <i>T</i><sub>2</sub><sup>*</sup>: ')
fwhm_w = winput( width=50, text='1.0', pos=scene.caption_anchor, 
                 bind=FWHM_select )

scene.append_to_caption(' <hr width="530" align="left">')

# Row 2 ----------------------------------------
scene.append_to_caption('  PROGRAM CONTROL - ')
button(text="Run", pos=scene.caption_anchor, bind=Run)
scene.append_to_caption(' ')
button(text="Reset M", pos=scene.caption_anchor, bind=Reset_M)
scene.append_to_caption(' ')
button(text="Clear time", pos=scene.caption_anchor, bind=Clear_plot)
scene.append_to_caption(' ')
frame_button = button(text="<b>Rotating</b>/Lab", pos=scene.caption_anchor,
                      bind=Rotate_frame)

# Row 3 -----------------------------------------
scene.append_to_caption('\n  Simulation speed:')
sl_sim_speed = slider(min=4, max=200, value=100, step=2, top=10, length=150,
                      right=15, bind=Set_sim_speed)
wt_sim_speed = wtext(text='{:.0f}%'.format(sl_sim_speed.value))
scene.append_to_caption('  ')
RF_frame_checkbox = checkbox(text="RF frame only", pos=scene.caption_anchor, 
                            checked=False, bind=Rotating_frame_select)           

scene.append_to_caption(' <hr width="530" align="left">')

# Row 4 -----------------------------------------
scene.append_to_caption('  CONTINUOUS RF: ')
button(text="On/<b>OFF</b>", pos=scene.caption_anchor, bind=RF_start)
scene.append_to_caption('   Align: ')
RF_menu = menu(choices=['x-axis','y-axis','oscillating'],
               pos=scene.caption_anchor, bind=RF_select )

scene.append_to_caption(' <hr width="530" align="left">')

# Row 5 -----------------------------------------
scene.append_to_caption('  PULSED RF: ')
button(text="Trigger", pos=scene.caption_anchor, bind=Sequence_start)
scene.append_to_caption('  ')
Instant_checkbox = checkbox(text="'Instant' RF", pos=scene.caption_anchor, 
                            checked=True, bind=Instant_select)               
scene.append_to_caption('  N<sub>B</sub>: ')
NumB_pulse_menu = menu( choices=['0','1','2','3','5','8','10','15','20'],
                        selected='0', pos=scene.caption_anchor, 
                        bind=NumB_pulse_select )
scene.append_to_caption(' Delay Time: ')
Interval_menu = menu( choices=['0.5','1','2','4','5','6','8','10','15','18',
                               '20','25','30','40','50'],
                      selected='5', pos=scene.caption_anchor,
                      bind=Interval_select )

# Row 6 -----------------------------------------
scene.append_to_caption('\n ')
Apulse_checkbox = checkbox(text='A pulse: ', pos=scene.caption_anchor,
                           checked=True, bind=Apulse_on)
Apulse_menu = menu( choices=['pi/2','pi','pi+error','optimal'], selected='pi/2',
                    pos=scene.caption_anchor, bind=Apulse_select )                   
scene.append_to_caption('    Align: ')
A_RF_menu = menu( choices=['x-axis','y-axis','oscillating'], selected='x-axis',
                  pos=scene.caption_anchor, bind=Apulse_select )

# Row 7 ----------------------------------------
scene.append_to_caption('\n ')
Bpulse_checkbox = checkbox(text='B pulse: ', pos=scene.caption_anchor,
                           checked=False, bind=Bpulse_on)
Bpulse_menu = menu( choices=['pi/2','pi','pi+error', 'optimal'], selected='pi',
                    pos=scene.caption_anchor, bind=Bpulse_select )
scene.append_to_caption('    Align: ')
B_RF_menu = menu( choices=['x-axis','y-axis','oscillating'], selected='x-axis', 
                  pos=scene.caption_anchor, bind=Bpulse_select )
scene.append_to_caption('  pi error(%): ')
w_delta = winput( width=30, text='20', pos=scene.caption_anchor, 
                 bind=delta_select)


scene.append_to_caption('<hr width="530" align="left">')


# Row 8/9 ---------------------------------------
scene.append_to_caption('  RF amplitude <i>&omega;</i><sub>1</sub>:\n')
sl_h1 = slider(min=0.01, max=1.2, value=0.4, length=250, bind=Set_h1, left=25,
               right=15)
wt_h1 = wtext(text='{:1.3f}'.format(sl_h1.value))
scene.append_to_caption(' <i>&omega;</i><sub>0</sub>')

# Row 10/11 -------------------------------------
scene.append_to_caption('\n  RF detuning &Delta;<i>&omega;</i>:\n')
sl_dW = slider(min=-2.5, max=2.5, value=0.0, step=0.01, length=250, 
               bind=Set_Delta_omega, left=25, right=15)
wt_dW = wtext(text='{:1.3f}'.format(sl_dW.value))
scene.append_to_caption(' <i>&omega;</i><sub>1</sub> (')
wt_dW_0 = wtext(text='{:1.3f}'.format(sl_dW.value))
scene.append_to_caption(' <i>&omega;</i><sub>0</sub>)')

# Row 12 ----------------------------------------
scene.append_to_caption('<hr width="530" align="left">')
scene.append_to_caption('  <i>T</i><sub>1</sub> (half-life): ')
T_menu = menu(choices=['4.0','8.0','10.0','20.0','50.0','100.0','200.0',
                       '500.0','1000.0','Infinite'],
              selected='Infinite', bind=T_select)

scene.append_to_caption('   <i>T</i><sub>2</sub> (% of <i>T</i><sub>1</sub>): ')
T2_menu = menu(choices=['100','90','50','10','5','1',], selected='100',
               bind=T2_select )

# Row 13 and plot area --------------------------
scene.append_to_caption('<hr width="530" align="left">')
scene.append_to_caption('  ')
button(text="Start plot", pos=scene.caption_anchor, bind=Stop_plot)
scene.append_to_caption('  ')
show_signal_checkbox = checkbox(bind=Show_signal, 
                                text='<font color="blue"><i>M</i><sub>&perp;</sub></font>  ',
                                checked=True)
show_Mz_checkbox = checkbox(bind=Show_Mz, 
                   text='<font color="gray"><i>M<sub>z</sub></i></font>  ',
                   checked=False)
show_Mr_ave_checkbox = checkbox(bind=Show_Mr_ave, 
                   text='<font color="orange">&langle;<i>M<sub>&perp;</sub></i>&rangle;</font>  ',
                   checked=False)
show_Mx_checkbox = checkbox(bind=Show_Mx,
                        text='<font color="red"><i>M<sub>x</sub></i></font>  ',
                        checked=False)
show_My_checkbox = checkbox(bind=Show_My, 
                        text='<font color="green"><i>M<sub>y</sub></i></font>',
                        checked=False)
scene.append_to_caption('  \n\n')

# -----------------------------------------------
# Initialize the scene
# -----------------------------------------------

# Add labels of time and time to next B pulse
wt_time = label(text='Time: {:.1f}'.format(t), pixel_pos=True,
                align='right', pos=vector(160,60,0))
wt_next_pulse = label(text=' ', pos=vector(180,60,0), align='left',
                      pixel_pos=True, box=False, opacity=0)

# Draw the axes on the canvas
rotating_frame = draw_axes()

# Place image objects of net magnetization and field
M = arrow(pos=vector(0,0,0), axis=m0, shaftwidth=0.3, round=True,
          color=color.blue)
B = arrow(pos=vector(0,0,0), axis=Om0/gamma, shaftwidth=0.3,
          round=True, color=color.green)

# Initializations here because most depend on defaults of controls.
#
# Initial values of loop variables
phi1 = (Domega()+omega0)*t # phase of RF
dphi1 = 0.0
pulse_time = 0.0 # Time below which RF is on
OMEGA = vector(0,0,omega0) # initial field vector
m = chi*OMEGA/gamma # initial magnetization
T_max = 1005.0/log(2) # sets value above which is "infinite"
T = T0/log(2) # Decay time constants 
T2 = T

# 
# Initial control states
#
running = False
in_rotating_frame = True
RF_frame_only = False
RF_on = False
RFtype='x-axis'
N_start = 0
sequence_state = 'Fire_A'
plotting = False
phi_latch = True
is_arrow = False
prog_ctl = 'init'
newspins = int(Spin_number.selected)
FWHM = sqrt(2.0*log(2.0))/float(fwhm_w.text)
Npulse = int(NumB_pulse_menu.selected)
Apulse_type =(Apulse_menu.selected, A_RF_menu.selected)
Bpulse_type =(Bpulse_menu.selected, B_RF_menu.selected)
sequence_interval = float(Interval_menu.selected)
t_nextPulse = 0.0
t_boundary = 0.0
instant_pulse = Instant_checkbox.checked
Apulse_yes = Apulse_checkbox.checked
Bpulse_yes = Bpulse_checkbox.checked
delta = 1.0 + float(w_delta.text)*0.01


# calc_OMEGA() is a central function, so it is declared here to better
# see how it affect program execution

def calc_OMEGA():
    '''
    Calculates the omega vector at any instant, depending on the current state
    of the system.  Also has a number of side effects for the control variables.
    Also calculates OMEGA_lab vector in lab frame for decay.
    '''
    global phi_latch, phi1, pulse_time, RFtype
    global Npulse, N_start, sequence_state, t_nextPulse
    global Apulse_type, Bpulse_type, tip
    global t_boundary, Tpulse

    this_OMEGA = vector(0.,0.,0.) # Local omega
    dt = current_dt # sets/resets in case dt adjusted below
    
    if Npulse > 0: 
        # If true, we are in a pulse sequence
                    
        if sequence_state == 'Fire_A':
            
            N_start = Npulse # Save initial pulse number
            
            if Apulse_yes: # Set pulse time and type to allow pulse
                pulse_time = t + Pulse_length(Apulse_type[0])
                RFtype = Apulse_type[1]

            Npulse = Npulse - 1
            if Npulse > 0: # Set time to next pulse, if it exists
                t_nextPulse = Tpulse[N_start-Npulse]
                t_boundary = t_nextPulse - dt*(1.0001)
                wt_next_pulse.text='B pulse {:d} at {:.1f}'.format(N_start-Npulse,t_nextPulse)
                # print(sequence_state, 'B in ', '{:.3f}'.format(t_nextPulse))
            else:
                wt_next_pulse.text=' '
            sequence_state = 'Idle'
            
        elif sequence_state == 'Fire_B':

            if Bpulse_yes: # Set time and type to allow next pulse
                pulse_time = t + Pulse_length(Bpulse_type[0])
                RFtype = Bpulse_type[1]
        
            Npulse = Npulse - 1
            if Npulse > 0: # Set time to next pulse, if i exists
                t_nextPulse = Tpulse[N_start-Npulse]
                t_boundary = t_nextPulse - dt*1.0001
                wt_next_pulse.text='B pulse {:d} at {:.1f}'.format(N_start-Npulse,t_nextPulse)
                # print(sequence_state, 'B in ', '{:.3f}'.format(t_nextPulse))
            else: # Otherwise, pulses are done.
                wt_next_pulse.text=' '
            sequence_state = 'Idle'

        elif sequence_state == 'Idle':
            # Due to the state-machine design of this function, we anticipate
            # the next pulse one step before it is due: t_nextPulse - dt
            if Npulse > 0 and t_plot > t_boundary: # small extra ensures boundary cross
                sequence_state = 'Fire_B'
            else:
                sequence_state = 'Idle'
    
    if RF_on or (t < pulse_time):
        # If t < pulse_time or RFon is true, add extra components to OMEGA 
        instant = instant_pulse  # instant prevents spins from evolving under local field

        if t < pulse_time and t+dt > pulse_time:
            # Shortens dt to account for step-size mismatch with pulse length
            dt = pulse_time - t
        
        if phi_latch:  # phi_latch resets the phase parameter to zero when using the
            phi1 = 0.0 # oscillating RF

        this_OMEGA.z = Domega()
        this_OMEGA.x, this_OMEGA.y = RF_function( phi1, RFtype )
        phi_latch = False
    else:
        # If there is no pulse, then the effective field is zero, and we are in the
        # resonant frame of reference.  The only evolution is for field inhomogeneity.
        if spin_0_trail: tip.clear_trail()
        #if t_plot > 0 and t_plot < t_nextPulse + default_dt and t_plot + 2*default_dt > t_nextPulse:
        #    print('t_plot = {:.5f},  t = {:.5f},  t_next = {:.5f}'.format(t_plot,t,t_nextPulse))
        
        instant = False # clear "instant" when there is no RF
        
        phi_latch = True  # resets the latch for the next pulse.
        
        this_OMEGA.x, this_OMEGA.y, this_OMEGA.z = 0.0, 0.0, 0.0
        if RF_frame_only:
            this_OMEGA.z = Domega()
 
    this_OMEGA_lab = vector(this_OMEGA.x, this_OMEGA.y, omega0)

    return dt, this_OMEGA, instant, this_OMEGA_lab

# ----------------------------------------------
# Main Loop
# ----------------------------------------------

# Initialize plot
gd, f1, mz, Mx, My, mr_ave = makeplot()

# Set up and use a trail on spin[0] when RF is on.  Helps visualise effect of 
# pulse on ensemble
spin_0_trail = False
if spin_0_trail:
    trail_points=300
    tip = sphere(radius=0.00, pos=m0, make_trail=True, color=color.blue,
                 visible=True, retain=trail_points)

while True:
    rate(ticks_per_s*pct_speed)
    
    if prog_ctl == 'init':
        Nspins = spin_initialize(Hwidth = FWHM, N = newspins)
        prog_ctl = 'operate'
        
    elif prog_ctl == 'operate':

        dt, OMEGA, instant, OMEGA_lab = calc_OMEGA()
        #if dt != default_dt:
        #    print('dt at {:.2f}%'.format(100*(dt)/default_dt))
        if is_arrow:
            for spin, em in zip(spins, lil_ems):
                em.axis = spin
        else:
            for spin, em in zip(spins, lil_ems):
                em.pos = spin
        ## Adjust the B (field) axis to effective value when in non-rotating frame
        if not in_rotating_frame:
            B.axis = OMEGA_lab/gamma
        else:
            B.axis = OMEGA/gamma
        m, m_ave_perp = ave_vector(spins)
        M.axis = m
        if spin_0_trail: tip.pos = lil_ems[0].pos
        
        if running:
            if not in_rotating_frame: # Update the camera position 
                scene.forward = scene.forward.rotate(angle=dphi1, axis=Om0)
            for i in range(len(spins)):
                mm = spins[i]
                if not instant: 
                    hh = OMEGA + omega_local[i]
                else: 
                    hh = OMEGA
                mm = mm.rotate(angle=-hh.mag*dt, axis=hh)                
                if T <= T_max and instant == False:
                     mm.z -= ((mm - (chi*(OMEGA_lab)/gamma)).z/T)*dt
                     mm.x -= ((mm - (chi*(OMEGA_lab)/gamma)).x/T2)*dt
                     mm.y -= ((mm - (chi*(OMEGA_lab)/gamma)).y/T2)*dt
                spins[i] = mm
            dphi1 = (-Domega()+omega0)*dt
            phi1 = phi1 + dphi1
            t = t + dt
            if not instant:
                t_plot = t_plot + dt
            #print('t {:.3f}, t_plot {:.3f}'.format(t,t_plot))
            wt_time.text='Time: {:.1f}'.format(t_plot)
            m_signal = sqrt(m.x**2 + m.y**2)
            if plotting: 
                f1.plot([t_plot , m_signal])
                mz.plot([t_plot, m.z])
                mr_ave.plot([t_plot, m_ave_perp])
                Mx.plot([t_plot, m.x])
                My.plot([t_plot, m.y])