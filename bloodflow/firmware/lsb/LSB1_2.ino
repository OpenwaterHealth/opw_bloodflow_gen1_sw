//Laser Safety Board
// By Peter Bradbury
// Rev 1

#define DIP 5
#define COIL 10
#define LED A5
#define PWR A2
#define BLED A3
#define HRTBT_OUT 11
#define HRTBT_IN A4
#define SR_FBK A0


int PD = 2;
int DIPs = 0;
int LEDstate = 0;
int flag1 = 0;
unsigned long duration;
unsigned long previousMillis = 0;
const long interval = 1000;
int PWRraw;
float PWRv;
int flag2 = 0;
int flag3 = 1;
unsigned long HRTBT_duration;
float SR_FBKv;

unsigned long timer = 0;
unsigned long endTimer = 0;

void setup() {
  pinMode(DIP, INPUT);
  pinMode(COIL, OUTPUT);
  pinMode(LED, OUTPUT);
  pinMode(PD, INPUT);
  pinMode(BLED, OUTPUT);
  pinMode(HRTBT_OUT, OUTPUT);
  pinMode(HRTBT_IN, INPUT);
  pinMode(SR_FBK, INPUT);
  //pinMode(PWR, INPUT);
  
  Serial1.begin(115200);
  delay(10);
  digitalWrite(COIL, HIGH);
//delay(10000);
//    while (!Serial)
//  {
//    ; // wait for serial port to connect. Needed for Leonardo only
//  }
  Serial.println("setup");
analogWrite(HRTBT_OUT, 125);
}

void loop() {

  duration = pulseIn(PD, LOW);
  Serial.print("Pulse Width: ");
  Serial.print(duration);
  Serial.println(" uS");
//  PWRraw = analogRead(PWR);
//  delay(50);
//  PWRraw = analogRead(PWR);
//  delay(50);
//  PWRv = PWRraw*(5.0/1023.0);
//   Serial.print("PWR voltage: ");
//  Serial.print(PWRv);
//  Serial.println(" V");
  
  unsigned long currentMillis = millis();

 if (currentMillis - previousMillis >= interval) {
    // save the last time you blinked the LED
    previousMillis = currentMillis;

    // if the LED is off turn it on and vice-versa:
    if (LEDstate == LOW) {
      LEDstate = HIGH;
    } else {
      LEDstate = LOW;
    }

    // set the LED with the ledState of the variable:
    digitalWrite(LED, LEDstate);
  
  }
  
 //DIPs = digitalRead(DIP);
//if (DIPs == HIGH){
//  digitalWrite (COIL, HIGH);
//  Serial.println("high");
//  
//}
//
//if (DIPs == LOW){
//  digitalWrite (COIL, LOW);
//  Serial.println("LOW");
//}
if (duration <= 430){
  digitalWrite (COIL, HIGH);
    //Serial.println("HIGH");
    Serial.print ("HIGH");
}


//if (duration <= 310 && duration >= 150){ 
//  digitalWrite (COIL, LOW);
//  Serial.println("LOW");
//}
if (duration >= 440){
   digitalWrite (COIL, LOW);
  Serial.println("LOW");
}

if (duration == 0){
  if (flag2 == 0){
  timer = millis();
  flag2 = 1;
  }
  endTimer = millis();
  if (endTimer - timer >= 200000){
    flag2 = 0;
    flag3 = 1;
  }
}

if (flag3 == 1){
  HRTBT_duration = pulseIn(HRTBT_IN, HIGH);
    SR_FBKv = analogRead(SR_FBK);
  flag3 = 0;
}

if (HRTBT_duration == 0){
  digitalWrite (COIL, LOW);
}

if (COIL == LOW && SR_FBKv >= 3.5){
  Serial.println("MCU SR control malfunciton or SR malfunction");
  Serial.println("MCU requesting SR coil low, but SR coil is high or contacts welded");
}

if (COIL == HIGH && SR_FBKv <= 1.5){
    Serial.println("MCU SR control malfunciton or SR malfunction");
  Serial.println("MCU requesting SR coil HIGH, but SR coil is LOW or contacts welded");
}

}
