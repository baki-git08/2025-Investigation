/******************************************************************************
 * Tailored Arduino Code for a Rotary Inverted Pendulum
 * * Control Scheme: Feedback Linearization with a PID Controller.
 * Hardware: Arduino Mega 2560, 1 Motor (L298N-style driver), 2 Encoders.
 * - Encoder 1 (Arm): Uses hardware interrupts (Pins 2, 3).
 * - Encoder 2 (Pendulum): Uses pin change interrupts (Pins 20, 21).
 * * Generated based on user-provided dynamics and parameters.
 ******************************************************************************/

// ----------------------- LIBRARIES ---------------------------------
#if defined(ESP32)
#define ISR_ATTR IRAM_ATTR
#else
#define ISR_ATTR
#endif

// ----------------------- USER SETTINGS ---------------------------
// --- Physical Parameters ---
const float g = 9.81;    // Gravity [m/s^2]
const float m1 = 0.25;    // Mass of the rotary arm
const float m2 = 0.032;    // Mass of the pendulum [kg]
const float L1 = 0.18;   // Length of the arm [m]
const float J1 = m1*(L1*L1);
const float L2 = 0.115;
const float L2 = l2/2;
const float J2 = m2/3*(L2*L2);  // Inertia of the pendulum [kg*m^2]
const float C1 = 0.01;    // Damping/friction coefficient for the arm

// --- STATE FEEDBACK CONTROLLER GAINS ---

// --- PID Controller Gains ---
const float Kp = 120.0;
const float Ki = 0.2;
const float Kd = 25.0;

// --- Control Target ---
float target_angle = 0.0; // The desired angle for the pendulum (0 = upright)

// --- Motor Pins (L298N-style driver) ---
const int motorPin1 = 8;  // Controls one side of the motor (e.g., IN1)
const int motorPin2 = 9;  // Controls the other side of the motor (e.g., IN2)
const int enablePin = 10; // PWM speed control pin (e.g., ENA)

// --- Encoder 1 Pins (Arm - uses Hardware Interrupts) ---
const uint8_t ARM_ENC_A = 20;
const uint8_t ARM_ENC_B = 21;

// --- Encoder 2 Pins (Pendulum - uses Pin Change Interrupts) ---
const uint8_t PENDULUM_ENC_A = 2;
const uint8_t PENDULUM_ENC_B = 3;

// --- Encoder Resolution ---
const long PPR = 2048; // Pulses Per Revolution (same for both encoders)

// --- Control System Limits ---
const float MAX_TORQUE = 5.0;     // Max torque for PWM mapping [N*m]
const float INTEGRAL_LIMIT = 100.0; // Anti-windup limit for the integral term

// ----------------------- GLOBAL VARIABLES ------------------------
// State variables (the "x" vector)
float armAngle_rad = 0.0;           // x1: Arm angle
float armVelocity_rad_s = 0.0;      // x2: Arm angular velocity
float pendulumAngle_rad = 0.0;      // x3: Pendulum angle
float pendulumVelocity_rad_s = 0.0; // x4: Pendulum angular velocity

// Raw encoder positions (updated by interrupts)
volatile long armEncoderPosition = 0;
volatile long pendulumEncoderPosition = 0;

// PID controller variables
float error_integral = 0.0;
unsigned long last_time = 0;

// Variables for speed calculation
long lastArmPosition = 0;
long lastPendulumPosition = 0;
unsigned long lastArmUpdateTime = 0;
unsigned long lastPendulumUpdateTime = 0;

// Variables for quadrature decoding in ISRs
volatile int lastArmEncoded = 0;
volatile int lastPendulumEncoded = 0;

// Define the absolute angle limits for the pendulum.
const float PENDULUM_ANGLE_LIMIT = (PI / 4.0); // +/- 45 degrees  

// Define the limits for torque
const float MAXIMUM_TORQUE = 5;
const float MINIMUM_TORQUE = -5;

// ----------------------- PHYSICS & CONTROL FUNCTIONS ----------------

/**
 * @brief Calculates the Lf2h and LgLfh terms from the system's dynamic model.
 * @param x2 Current arm angular velocity.
 * @param x3 Current pendulum angle.
 * @param x4 Current pendulum angular velocity.
 * @param Lf2h Reference to store the calculated Lf2h.
 * @param LgLfh Reference to store the calculated LgLfh.
 **/
void calculateParameters(float &x2, float &x3, float &x4, float &Lf2h, float &LgLfh) {
    float cos_x3 = cos(x3);
    float sin_x3 = sin(x3);

    // Common denominator for both equations to save computation
    float denominator = (L1*L1 * l2*l2 * m2*m2 * cos_x3*cos_x3) + 
                        (L1*L1 * l2*l2 * m2*m2) + 
                        (J2 * L1*L1 * m2) + 
                        (J1 * l2*l2 * m2) + // Corrected from your notes: J1*l2^2*m2
                        (J1 * J2);

    // Calculate Lf2h
    float lf2h_numerator = l2*m2 * (l2*m2*cos_x3*sin_x3*L1*L1*x4*x4 + 
                                    g*m2*sin_x3*L1*L1 + 
                                    C1*x2*cos_x3*L1 + 
                                    J1*g*sin_x3);
    Lf2h = lf2h_numerator / denominator;

    // Calculate LgLfh
    float lglfh_numerator = -L1 * l2 * m2 * cos_x3;
    LgLfh = lglfh_numerator / denominator;
}

/**
 * @brief Controls the motor based on a calculated torque value.
 * @param torque The desired output torque in N*m.
 */
void controlMotor(float torque) {

    int pwmValue = map(abs(torque), 0, MAX_TORQUE, 0, 255);

    if (pwmValue > 255) {
        pwmValue = 255;
    }

    // Set motor direction based on the sign of the torque
    if (torque < 0) {
        digitalWrite(motorPin1, HIGH);
        digitalWrite(motorPin2, LOW);
    } else {
        digitalWrite(motorPin1, LOW);
        digitalWrite(motorPin2, HIGH);
    }
    analogWrite(enablePin, pwmValue);
}

// ----------------------- HELPER & SENSOR FUNCTIONS ---------------

/**
 * @brief Converts raw encoder position to radians.
 */
float positionToRadians(long position) {
    return (float)position * (2.0 * PI) / (float)PPR;
}

/**
 * @brief Calculates arm angular velocity in rad/s.
 */
void calculateArmSpeed() {
    unsigned long currentTime = micros();
    unsigned long timeDelta = currentTime - lastArmUpdateTime;

    if (timeDelta > 0) {
        long positionDelta = armEncoderPosition - lastArmPosition;
        armVelocity_rad_s = ((float)positionDelta / (float)PPR) * (2.0 * PI) / ((float)timeDelta / 1000000.0);
        lastArmPosition = armEncoderPosition;
        lastArmUpdateTime = currentTime;
    }
}

/**
 * @brief Calculates pendulum angular velocity in rad/s.
 */
void calculatePendulumSpeed() {
    unsigned long currentTime = micros();
    unsigned long timeDelta = currentTime - lastPendulumUpdateTime;

    if (timeDelta > 0) {
        long positionDelta = pendulumEncoderPosition - lastPendulumPosition;
        pendulumVelocity_rad_s = ((float)positionDelta / (float)PPR) * (2.0 * PI) / ((float)timeDelta / 1000000.0);
        lastPendulumPosition = pendulumEncoderPosition;
        lastPendulumUpdateTime = currentTime;
    }
}

// ----------------------- INTERRUPT SERVICE ROUTINES --------------

/**
 * @brief ISR for the Arm Encoder (Hardware Interrupt).
 */
void updateArmEncoder() {
    int MSB = digitalRead(ARM_ENC_A);
    int LSB = digitalRead(ARM_ENC_B);
    int encoded = (MSB << 1) | LSB;
    int sum = (lastArmEncoded << 2) | encoded;

    if (sum == 0b1101 || sum == 0b0100 || sum == 0b0010 || sum == 0b1011) armEncoderPosition++;
    if (sum == 0b1110 || sum == 0b0111 || sum == 0b0001 || sum == 0b1000) armEncoderPosition--;
    
    lastArmEncoded = encoded;
}

/**
 * @brief ISR for the Pendulum Encoder (Pin Change Interrupt).
 */
void updatePendulumEncoder() {
    int MSB = digitalRead(PENDULUM_ENC_A);
    int LSB = digitalRead(PENDULUM_ENC_B);
    int encoded = (MSB << 1) | LSB;
    int sum = (lastPendulumEncoded << 2) | encoded;

    if (sum == 0b1101 || sum == 0b0100 || sum == 0b0010 || sum == 0b1011) pendulumEncoderPosition++;
    if (sum == 0b1110 || sum == 0b0111 || sum == 0b0001 || sum == 0b1000) pendulumEncoderPosition--;
    
    lastPendulumEncoded = encoded;
}

// ----------------------- MAIN PROGRAM ----------------------------

void setup() {
    Serial.begin(115200);

    // Motor pin setup
    pinMode(motorPin1, OUTPUT);
    pinMode(motorPin2, OUTPUT);
    pinMode(enablePin, OUTPUT);

    // Encoder pin setup
    pinMode(ARM_ENC_A, INPUT_PULLUP);
    pinMode(ARM_ENC_B, INPUT_PULLUP);
    pinMode(PENDULUM_ENC_A, INPUT_PULLUP);
    pinMode(PENDULUM_ENC_B, INPUT_PULLUP);

    // Attach Arm Encoder to Hardware Interrupts
    attachInterrupt(digitalPinToInterrupt(ARM_ENC_A), updateArmEncoder, CHANGE);
    attachInterrupt(digitalPinToInterrupt(ARM_ENC_B), updateArmEncoder, CHANGE);
    attachInterrupt(digitalPinToInterrupt(PENDULUM_ENC_A), updatePendulumEncoder, CHANGE);
    attachInterrupt(digitalPinToInterrupt(PENDULUM_ENC_B), updatePendulumEncoder, CHANGE);

    // Serial.println("Rotary Inverted Pendulum Initialized.");
    last_time = micros();
}

void loop() {
    // --- 1. SENSE: Get the current state of the system ---
    armAngle_rad = positionToRadians(armEncoderPosition);
    pendulumAngle_rad = positionToRadians(pendulumEncoderPosition);
    calculateArmSpeed();
    calculatePendulumSpeed();

    // Time difference for PID integral and derivative
    unsigned long now = micros();
    float dt = (float)(now - last_time) / 1000000.0;
    last_time = now;

    // --- 2. CONTROL: Calculate the required torque ---
    float error = target_angle - pendulumAngle_rad;
    
    // Integral term with anti-windup
    error_integral += error * dt;
    if (error_integral > INTEGRAL_LIMIT) error_integral = INTEGRAL_LIMIT;
    if (error_integral < -INTEGRAL_LIMIT) error_integral = -INTEGRAL_LIMIT;

    // Derivative term (using velocity is more stable than change in error)
    float error_dot = -pendulumVelocity_rad_s;

    float error = target_angle - pendulumAngle_rad;
    error_integral += error * dt_s;
    if (error_integral > INTEGRAL_LIMIT) error_integral = INTEGRAL_LIMIT;
    if (error_integral < -INTEGRAL_LIMIT) error_integral = -INTEGRAL_LIMIT;
    float error_dot = -pendulumVelocity_rad_s; // Derivative term on measurement

    // PID CONTROLLER
    float u_new = (Kp * error) + (Ki * error_integral) + (Kd * error_dot);

    // b) Feedback Linearization calculates the final torque
    float Lf2h, LgLfh, torque;
    calculateParameters(armVelocity_rad_s, pendulumAngle_rad, pendulumVelocity_rad_s, Lf2h, LgLfh);

    // Avoid division by zero if LgLfh is ever zero (e.g., if pendulum is horizontal)
    if (abs(LgLfh) > 1e-6) {
        torque = (u_new - Lf2h) / LgLfh;
    }

    // --- 3. ACTUATE: Apply the calculated torque to the motor ---
    // Check if the pendulum's absolute angle is within the safe operating zone.
    if (abs(pendulumAngle_rad) <= PENDULUM_ANGLE_LIMIT) {
        
        // Clamp the torque to its maximum safe value.
        if (torque > MAXIMUM_TORQUE) {
            torque = MAXIMUM_TORQUE;
        } else if (torque < MINIMUM_TORQUE) {
            torque = MINIMUM_TORQUE;
        }
        // Send the final, clamped torque to the motor.
        controlMotor(torque);
    } else {
        controlMotor(0.0);
    }

    /*  // --- 4. DEBUG: Print key values to Serial Plotter/Monitor ---
    Serial.print("\tArmA(rad):"); Serial.print(armAngle_rad, 2);
    Serial.print("\tPenA(rad):"); Serial.print(pendulumAngle_rad, 2);
    Serial.print("\tTorque:"); Serial.print(torque, 2);*/
}

