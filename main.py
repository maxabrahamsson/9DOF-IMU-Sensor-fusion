import time,math
from i2clibraries import i2c_adxl345,i2c_hmc5883l,i2c_itg3205
from array import array
# CALIBRATION VALUEs
ACCEL_X_MIN = -289
ACCEL_X_MAX = 294
ACCEL_Y_MIN = -268
ACCEL_Y_MAX = 288
ACCEL_Z_MIN = -294
ACCEL_Z_MAX = 269

MAGN_X_MIN = -600
MAGN_X_MAX = 600
MAGN_Y_MIN = -600
MAGN_Y_MAX = 600
MAGN_Z_MIN = -600
MAGN_Z_MAX = 600

GYRO_AVERAGE_OFFSET_X = 23.85
GYRO_AVERAGE_OFFSET_Y = -53.41
GYRO_AVERAGE_OFFSET_Z = -15.32

CALIBRATION__MAGN_USE_EXTENDED = True
magn_ellipsoid_center = [3.80526, -16.4455, 87.4052]
magn_ellipsoid_transform = [[0.970991, 0.00583310, -0.00265756], [0.00583310, 0.952958, 2.76726e-05], [-0.00265756, 2.76726e-05, 0.999751]]
##
GRAVITY = 256.0 # "1G reference" used for DCM filter and accelerometer calibration

ACCEL_X_OFFSET = ((ACCEL_X_MIN + ACCEL_X_MAX) / 2.0)
ACCEL_Y_OFFSET = ((ACCEL_Y_MIN + ACCEL_Y_MAX) / 2.0)
ACCEL_Z_OFFSET = ((ACCEL_Z_MIN + ACCEL_Z_MAX) / 2.0)
ACCEL_X_SCALE = (GRAVITY / (ACCEL_X_MAX - ACCEL_X_OFFSET))
ACCEL_Y_SCALE = (GRAVITY / (ACCEL_Y_MAX - ACCEL_Y_OFFSET))
ACCEL_Z_SCALE = (GRAVITY / (ACCEL_Z_MAX - ACCEL_Z_OFFSET))

MAGN_X_OFFSET = ((MAGN_X_MIN + MAGN_X_MAX) / 2.0)
MAGN_Y_OFFSET = ((MAGN_Y_MIN + MAGN_Y_MAX) / 2.0)
MAGN_Z_OFFSET = ((MAGN_Z_MIN + MAGN_Z_MAX) / 2.0)
MAGN_X_SCALE = (100.0 / (MAGN_X_MAX - MAGN_X_OFFSET))
MAGN_Y_SCALE = (100.0 / (MAGN_Y_MAX - MAGN_Y_OFFSET))
MAGN_Z_SCALE = (100.0 / (MAGN_Z_MAX - MAGN_Z_OFFSET))

def TO_RAD(x):
  return (x * 0.01745329252)

def GYRO_SCALED_RAD(x): 
  return (x * TO_RAD(GYRO_GAIN)) 
GYRO_GAIN = 0.06957 # Same gain on all axes

Kp_ROLLPITCH = 0.02
Ki_ROLLPITCH = 0.00002
Kp_YAW = 1.2
Ki_YAW = 0.00002


def Normalize():
  error=0
  temporary=[[0,0,0],[0,0,0],[0,0,0]]
  renorm=0
  
  error= -Vector_Dot_Product(DCM_Matrix[0],DCM_Matrix[1])*.5 #eq.19

  temporary[0] = Vector_Scale(DCM_Matrix[1], error) #eq.19
  temporary[1] = Vector_Scale(DCM_Matrix[0], error) #eq.19
  
  temporary[0]=Vector_Add(temporary[0], DCM_Matrix[0])#eq.19
  temporary[1]=Vector_Add(temporary[1], DCM_Matrix[1])#eq.19
  
  temporary[2]=Vector_Cross_Product(temporary[0],temporary[1]) # c= a x b #eq.20
  
  renorm= 0.5 *(3 - Vector_Dot_Product(temporary[0],temporary[0])) #eq.21
  DCM_Matrix[0] = Vector_Scale(temporary[0], renorm)
  
  renorm= 0.5 *(3 - Vector_Dot_Product(temporary[1],temporary[1])) #eq.21
  DCM_Matrix[1] = Vector_Scale(temporary[1], renorm)
  
  renorm= 0.5 *(3 - Vector_Dot_Product(temporary[2],temporary[2])) #eq.21
  DCM_Matrix[2] = Vector_Scale(temporary[2], renorm)

def Drift_correction():
  global Omega_I,Accel_Vector,GRAVITY,DCM_Matrix,Kp_ROLLPITCH,Ki_ROLLPITCH
  global MAG_Heading,Kp_YAW
  mag_heading_x=0
  mag_heading_y=0
  #Compensation the Roll, Pitch and Yaw drift. 
  Scaled_Omega_P=[0,0,0]
  Scaled_Omega_I=[0,0,0]
  Accel_magnitude=0
  Accel_weight=0
  
  
  #*****Roll and Pitch***************

  # Calculate the magnitude of the accelerometer vector
  Accel_magnitude = math.sqrt(Accel_Vector[0]*Accel_Vector[0] + Accel_Vector[1]*Accel_Vector[1] + Accel_Vector[2]*Accel_Vector[2]);
  Accel_magnitude = Accel_magnitude / GRAVITY; # Scale to gravity.
  # Dynamic weighting of accelerometer info (reliability filter)
  # Weight for accelerometer info (<0.5G = 0.0, 1G = 1.0 , >1.5G = 0.0)
  Accel_weight = min(max(1 - 2*abs(1 - Accel_magnitude),1),0)  #  

  errorRollPitch=Vector_Cross_Product(Accel_Vector,DCM_Matrix[2]) #adjust the ground of reference
  Omega_P=Vector_Scale(errorRollPitch,Kp_ROLLPITCH*Accel_weight)
  
  Scaled_Omega_I=Vector_Scale(errorRollPitch,Ki_ROLLPITCH*Accel_weight)
  Omega_I=Vector_Add(Omega_I,Scaled_Omega_I)
  
  #*****YAW***************
  # We make the gyro YAW drift correction based on compass magnetic heading
 
  mag_heading_x = math.cos(MAG_Heading)
  mag_heading_y = math.sin(MAG_Heading)
  errorCourse=(DCM_Matrix[0][0]*mag_heading_y) - (DCM_Matrix[1][0]*mag_heading_x)  #Calculating YAW error
  errorYaw=Vector_Scale(DCM_Matrix[2],errorCourse) #Applys the yaw correction to the XYZ rotation of the aircraft, depeding the position.
  
  Scaled_Omega_P=Vector_Scale(errorYaw,Kp_YAW) #.01proportional of YAW.
  Omega_P=Vector_Add(Omega_P,Scaled_Omega_P) #Adding  Proportional.
  
  Scaled_Omega_I=Vector_Scale(errorYaw,Ki_YAW)#.00001Integrator
  Omega_I=Vector_Add(Omega_I,Scaled_Omega_I)#adding integrator to the Omega_I


def Matrix_update():
  Gyro_Vector[0]=GYRO_SCALED_RAD(gyro[0]) #gyro x roll
  Gyro_Vector[1]=GYRO_SCALED_RAD(gyro[1]) #gyro y pitch
  Gyro_Vector[2]=GYRO_SCALED_RAD(gyro[2]) #gyro z yaw
  
  Accel_Vector[0]=accel[0]
  Accel_Vector[1]=accel[1]
  Accel_Vector[2]=accel[2]
    
  Omega=Vector_Add(Gyro_Vector, Omega_I)  #adding proportional term
  Omega_Vector=Vector_Add(Omega, Omega_P) #adding Integrator term
  

# Use drift correction
  Update_Matrix[0][0]=0
  Update_Matrix[0][1]=-G_Dt*Omega_Vector[2]#-z
  Update_Matrix[0][2]=G_Dt*Omega_Vector[1]#y
  Update_Matrix[1][0]=G_Dt*Omega_Vector[2]#z
  Update_Matrix[1][1]=0
  Update_Matrix[1][2]=-G_Dt*Omega_Vector[0]#-x
  Update_Matrix[2][0]=-G_Dt*Omega_Vector[1]#-y
  Update_Matrix[2][1]=G_Dt*Omega_Vector[0]#x
  Update_Matrix[2][2]=0

  Matrix_Multiply(DCM_Matrix,Update_Matrix,Temporary_Matrix) #a*b=c

  for x in range(0,3):
    for y in range(0,3):
      DCM_Matrix[x][y]+=Temporary_Matrix[x][y]

def Euler_angles():
  global pitch,roll,yaw
  pitch = -math.asin(DCM_Matrix[2][0])
  roll = math.atan2(DCM_Matrix[2][1],DCM_Matrix[2][2])
  yaw = math.atan2(DCM_Matrix[1][0],DCM_Matrix[0][0])

def Compass_Heading():
  cos_roll = math.cos(roll)
  sin_roll = math.sin(roll)
  cos_pitch = math.cos(pitch)
  sin_pitch = math.sin(pitch)
  
  # Tilt compensated magnetic field X
  mag_x = magnetom[0] * cos_pitch + magnetom[1] * sin_roll * sin_pitch + magnetom[2] * cos_roll * sin_pitch
  # Tilt compensated magnetic field Y
  mag_y = magnetom[1] * cos_roll - magnetom[2] * sin_roll
  # Magnetic Heading
  MAG_Heading = math.atan2(-mag_y, mag_x)

# Computes the dot product of two vectors
def Vector_Dot_Product(v1, v2):
  result = 0.0
  for c in range(0,3):
    result += v1[c] * v2[c]
  return result

# Computes the cross product of two vectors
# out has to different from v1 and v2 (no in-place)!
def Vector_Cross_Product(v1, v2):
  out=[0,0,0]
  out[0] = (v1[1] * v2[2]) - (v1[2] * v2[1])
  out[1] = (v1[2] * v2[0]) - (v1[0] * v2[2])
  out[2] = (v1[0] * v2[1]) - (v1[1] * v2[0])
  return out

# Multiply the vector by a scalar
def Vector_Scale(v, scale):
  out=[0,0,0]
  for c in range(0,3):
    out[c] = v[c] * scale
  return out

# Adds two vectors
def Vector_Add(v1, v2):
  out=[0,0,0]
  for c in range(0,3):
    out[c] = v1[c] + v2[c]
  return out

# Multiply two 3x3 matrices: out = a * b
# out has to different from a and b (no in-place)!
def Matrix_Multiply(a, b,out):
  for x in range(0,3): 
    for y in range(0,3):
      out[x][y] = a[x][0] * b[0][y] + a[x][1] * b[1][y] + a[x][2] * b[2][y]

# Multiply 3x3 matrix with vector: out = a * b
# out has to different from b (no in-place)!
def Matrix_Vector_Multiply(a, b,out):
  for x in range(0,3):
    out[x] = a[x][0] * b[0] + a[x][1] * b[1] + a[x][2] * b[2]

# Init rotation matrix using euler angles
def init_rotation_matrix(m, yaw, pitch, roll):
  c1 = math.cos(roll)
  s1 = math.sin(roll)
  c2 = math.cos(pitch)
  s2 = math.sin(pitch);
  c3 = math.cos(yaw)
  s3 = math.sin(yaw)

  m[0][0] = c2 * c3;
  m[0][1] = c3 * s1 * s2 - c1 * s3
  m[0][2] = s1 * s3 + c1 * c3 * s2

  m[1][0] = c2 * s3;
  m[1][1] = c1 * c3 + s1 * s2 * s3
  m[1][2] = c1 * s2 * s3 - c3 * s1

  m[2][0] = -s2
  m[2][1] = c2 * s1
  m[2][2] = c1 * c2

def TO_DEG(x):
  return (x * 57.2957795131)  # *180/pi
def output_angles():
  global yaw,pitch,roll
  print("#YPR= "+str(int(TO_DEG(yaw)))+" "+str(int(TO_DEG(pitch)))+" "+str(int(TO_DEG(roll))))

# Sensor variables
accel=[0,0,0]  # Actually stores the NEGATED acceleration (equals gravity, if board not moving).
accel_min=[0,0,0]
accel_max=[0,0,0]

magnetom=[0,0,0]
magnetom_min=[0,0,0]
magnetom_max=[0,0,0]
magnetom_tmp=[0,0,0]

gyro=[0,0,0]
gyro_average=[0,0,0]
gyro_num_samples = 0;

# DCM variables
MAG_Heading=0.0;
Accel_Vector= [0, 0, 0] # Store the acceleration in a vector
Gyro_Vector= [0,0,0] # Store the gyros turn rate in a vector
Omega_Vector = [0,0,0] # Corrected Gyro_Vector data
Omega_P= [0,0,0] # Omega Proportional correction
Omega_I= [0,0,0] # Omega Integrator
Omega = [0,0,0]
errorRollPitch = [0,0,0]
errorYaw = [0,0,0]
DCM_Matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
Update_Matrix = [[0, 1, 2],[3, 4, 5], [6, 7, 8]]
Temporary_Matrix = [[0, 0, 0],[0, 0, 0], [0, 0, 0]]

# Euler angles
yaw=0.0
pitch=0.0
roll=0.0

# DCM timing in the main loop\
timestamp=0
timestamp_old=0
G_Dt=0.0

# More output-state variables
curr_calibration_sensor = 0;
reset_calibration_session_flag = True;
num_accel_errors = 0;
num_magn_errors = 0;
num_gyro_errors = 0;



def read_sensors():
  gyro=itg3205.getAxes() # Read gyroscope
  accel=adxl345.getAxes() # Read accelerometer
  magnetom=hmc5883l.getAxes() # Read magnetometer
  #gyro=(0.2,0.1,0.2) # Read gyroscope
  #accel=(0.2,0.1,0.2) # Read accelerometer
  #magnetom=(0.2,0.1,0.2) # Read magnetometer
  print((gyro,accel,magnetom))
# Read every sensor and record a time stamp
# Init DCM with unfiltered orientation
# TODO re-init global vars?

def millis():
  return round(time.time() * 1000)

def reset_sensor_fusion():
  global DCM_Matrix,accel,timestamp,yaw,pitch,roll,MAG_Heading
  temp1=[0,0,0]
  temp2=[0,0,0]
  xAxis = [1.0, 0.0, 0.0]

  read_sensors()
  timestamp = millis()
  
  # GET PITCH
  # Using y-z-plane-component/x-component of gravity vector
  pitch = -math.atan2(accel[0], math.sqrt(accel[1] * accel[1] + accel[2] * accel[2]))
	
  # GET ROLL
  # Compensate pitch of gravity vector 
  temp1=Vector_Cross_Product( accel, xAxis)
  temp2=Vector_Cross_Product( xAxis, temp1)
  # Normally using x-z-plane-component/y-component of compensated gravity vector
  # roll = atan2(temp2[1], sqrt(temp2[0] * temp2[0] + temp2[2] * temp2[2]));
  # Since we compensated for pitch, x-z-plane-component equals z-component:
  roll = math.atan2(temp2[1], temp2[2])
  
  # GET YAW
  Compass_Heading()
  yaw = MAG_Heading
  
  # Init rotation matrix
  init_rotation_matrix(DCM_Matrix, yaw, pitch, roll)

# Apply calibration to raw sensor readings
def compensate_sensor_errors():
    accel[0] = (accel[0] - ACCEL_X_OFFSET) * ACCEL_X_SCALE
    accel[1] = (accel[1] - ACCEL_Y_OFFSET) * ACCEL_Y_SCALE
    accel[2] = (accel[2] - ACCEL_Z_OFFSET) * ACCEL_Z_SCALE

    # Compensate magnetometer error
#if CALIBRATION__MAGN_USE_EXTENDED == true
    for i in range(0,3):
      magnetom_tmp[i] = magnetom[i] - magn_ellipsoid_center[i]
    Matrix_Vector_Multiply(magn_ellipsoid_transform, magnetom_tmp, magnetom)
#else
    magnetom[0] = (magnetom[0] - MAGN_X_OFFSET) * MAGN_X_SCALE
    magnetom[1] = (magnetom[1] - MAGN_Y_OFFSET) * MAGN_Y_SCALE
    magnetom[2] = (magnetom[2] - MAGN_Z_OFFSET) * MAGN_Z_SCALE
#endif

    # Compensate gyroscope error
    gyro[0] -= GYRO_AVERAGE_OFFSET_X
    gyro[1] -= GYRO_AVERAGE_OFFSET_Y
    gyro[2] -= GYRO_AVERAGE_OFFSET_Z

# Reset calibration session if reset_calibration_session_flag is set
def check_reset_calibration_session():
  # Raw sensor values have to be read already, but no error compensation applied
  # Reset this calibration session?
  global reset_calibration_session_flag
  global accel_min,accel_max,accel
  global magnetom_max,magnetom_min,magnetom
  global gyro_average,gyro_num_samples
  if (reset_calibration_session_flag==False):
    return
  
  # Reset acc and mag calibration variables
  for i in range(0,3):
    accel_min[i] = accel_max[i] = accel[i]
    magnetom_min[i] = magnetom_max[i] = magnetom[i]

  # Reset gyro calibration variables
  gyro_num_samples = 0  # Reset gyro calibration averaging
  gyro_average[0] = gyro_average[1] = gyro_average[2] = 0.0
  
  reset_calibration_session_flag = False


def setup():
  #time.sleep(5)
  global adxl345,hmc5883l,itg3205

  adxl345 = i2c_adxl345.i2c_adxl345(3)
  hmc5883l = i2c_hmc5883l.i2c_hmc5883l(3)
  itg3205 = i2c_itg3205.i2c_itg3205(3)
  
  # Read sensors, init DCM algorithm
  time.sleep(2)  # Give sensors enough time to collect data
  reset_sensor_fusion()

# Main loop
def loop():
  global timestamp,timestamp_old,G_Dt
  # Time to read the sensors again?
  if((millis() - timestamp) >= 1000):
    timestamp_old = timestamp
    timestamp = millis()
    if (timestamp > timestamp_old):
      G_Dt = (timestamp - timestamp_old) / 1000.0
    else:
      G_Dt = 0

    # Update sensor readings
    read_sensors()
    # Apply sensor calibration
    compensate_sensor_errors()
  
    # Run DCM algorithm
    Compass_Heading() # Calculate magnetic heading
    Matrix_update()
    Normalize()
    Drift_correction()
    Euler_angles()
    
    output_angles()
    
    #print("loop time (ms) = ")
    #print(millis() - timestamp)

setup()
while True:
  loop()
