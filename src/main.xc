// COMS20001 - Cellular Automaton Farm - Initial Code Skeleton
// (using the XMOS i2c accelerometer demo code)

#include <platform.h>
#include <xs1.h>
#include <stdio.h>
#include "pgmIO.h"
#include "i2c.h"

#define IMHT 16                  //image height
#define IMWD 16                  //image width
#define WORKERS 8

typedef unsigned char uchar;      //using uchar as shorthand

port p_scl = XS1_PORT_1E;         //interface ports to orientation
port p_sda = XS1_PORT_1F;

#define FXOS8700EQ_I2C_ADDR 0x1E  //register addresses for orientation
#define FXOS8700EQ_XYZ_DATA_CFG_REG 0x0E
#define FXOS8700EQ_CTRL_REG_1 0x2A
#define FXOS8700EQ_DR_STATUS 0x0
#define FXOS8700EQ_OUT_X_MSB 0x1
#define FXOS8700EQ_OUT_X_LSB 0x2
#define FXOS8700EQ_OUT_Y_MSB 0x3
#define FXOS8700EQ_OUT_Y_LSB 0x4
#define FXOS8700EQ_OUT_Z_MSB 0x5
#define FXOS8700EQ_OUT_Z_LSB 0x6

void keepTime() {
    timer tmr;
    unsigned long t;
    tmr :> t;
    t /= 100000;
    //printf("%lu\n", t);
    unsigned int period = 100000000;
    select {
        case tmr when timerafter(period) :> void:
            t = ((period/100000) - t);
            //printf("%lu\n", t);
            break;
    }
    while (1) {
        timer tmr;
        select {
                case tmr when timerafter(period) :> void:
                    t += (period/100000);
                    //printf("%lu\n", t);
                    break;
            }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////
//
// Read Image from PGM file from path infname[] to channel c_out
//
/////////////////////////////////////////////////////////////////////////////////////////
void DataInStream(char infname[], chanend c_out)
{
  int res;
  uchar line[ IMWD ];
  printf( "DataInStream: Start...\n" );

  //Open PGM file
  res = _openinpgm( infname, IMWD, IMHT );
  if( res ) {
    printf( "DataInStream: Error openening %s\n.", infname );
    return;
  }

  //Read image line-by-line and send byte by byte to channel c_out
  for( int y = 0; y < IMHT; y++ ) {
    _readinline( line, IMWD );
    for( int x = 0; x < IMWD; x++ ) {
      c_out <: line[ x ];
    }
  }
  //Close PGM image file
  _closeinpgm();
  printf( "DataInStream: Done...\n" );
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////
//
// Start your implementation by changing this function to implement the game of life
// by farming out parts of the image to worker threads who implement it...
// Currently the function just inverts the image
//
/////////////////////////////////////////////////////////////////////////////////////////

void worker(chanend c) {
    //uchar slice[(IMHT/8) + 2][IMWD];
}

void slice(uchar result[IMHT/WORKERS + 2][IMWD/8], uchar image[IMHT][IMWD/8], int i) {
    if (i == 0){
        memcpy(image[IMHT-1], result[0], sizeof(image[IMHT-1]));
        int k = 1;
        for(int y = i; y <= (IMHT/WORKERS); y++){
            memcpy(image[y], result[k], sizeof(image[y]));
            k++;
        }
    } else if (i == (WORKERS - 1)){
        int k = 0;
        for(int y = i-1; y < (IMHT/WORKERS); y++){
            memcpy(image[y], result[k], sizeof(image[y]));
            k++;
        }
        memcpy(image[0], result[(IMWD/8)-1], sizeof(image[0]));
    } else {
        int k = 0;
        for(int y = i-1; y <= (IMHT/WORKERS); y++){
            memcpy(image[y], result[k], sizeof(image[y]));
            k++;
        }
    }
}

uchar getCell(uchar image[ IMHT ][ IMWD / 8 ], uchar x, uchar y){
    return image[y][x/8]&(1<<(x%8)) ? 1 : 0;
}

void setCell(uchar image[ IMHT ][ IMWD / 8], uchar x, uchar y, uchar v){
    if (v){
        image[y][x/8] |= (1<<(x%8));
    } else {
        (image[y][x/8]) &= ~(((uchar)1)<<(x%8));
    }
}

uchar checkCell(uchar image[ IMHT ][ IMWD / 8 ], uchar x, uchar y)
{
    uchar living = getCell(image, x, y);
    uchar neighbours = 0;
    uchar result = 0;
    uchar ym = (y-1); ym %= IMHT;
    uchar yp = (y+1); yp %= IMHT;
    uchar xm = (x-1); xm %= IMWD;
    uchar xp = (x+1); xp %= IMWD;
    //printf("ym:%u yp:%u xm:%u xp:%d\n", ym, yp, xm, xp);

    neighbours += getCell(image, xm, ym);
    neighbours += getCell(image, x, ym);
    neighbours += getCell(image, xp, ym);

    neighbours += getCell(image, xm, y);
    neighbours += getCell(image, xp, y);

    neighbours += getCell(image, xm, yp);
    neighbours += getCell(image, x, yp);
    neighbours += getCell(image, xp, yp);

    //printf("alive:%u n=%u\n", living, neighbours);
    //printf("l=%u", living);
    if(living != 0){
        result = (neighbours == 2 || neighbours == 3) ? 255 : 0;
    } else {
        result = (neighbours == 3) ? 255 : 0;
    }

    return result;
}

uchar ** getRow(uchar slice[IMHT/WORKERS + 2][IMWD/8], int iterations){
    uchar row[IMHT/WORKERS][IMWD/8];

    for (int i = 0; i < iterations; i++)
    {
      for( int y = 1; y < IMHT/WORKERS; y++ ) {   //go through all lines
        for( int x = 0; x < IMWD; x++ ) {
                setCell(imageNext, x, y, checkCell(imageCurrent, x, y));
        }
      }
      memcpy(imageCurrent, imageNext, sizeof(imageCurrent));
      printWorld(imageCurrent);
      printf("\n");
    }
}

void printWorld(uchar image[ IMHT ][ IMWD / 8]) {
    for( int y = 0; y < IMHT; y++ ) {
        for( int x = 0; x < IMWD; x++ ) {
          char p = (getCell(image, x, y)) ? 'X' : '.';
          printf( "%c ",p); //show image values
        }
        printf( "\n" );
      }
}

void generate(int iterations, uchar imageCurrent[IMHT][IMWD/8], uchar imageNext[IMHT][IMWD/8]){
    timer tmr;
    int startTime;
    tmr :> startTime;

    for (int i = 0; i < iterations; i++)
    {
      for( int y = 0; y < IMHT; y++ ) {   //go through all lines
        for( int x = 0; x < IMWD; x++ ) {
                setCell(imageNext, x, y, checkCell(imageCurrent, x, y));
        }
      }
      memcpy(imageCurrent, imageNext, sizeof(imageCurrent));
      printWorld(imageCurrent);
      printf("\n");
    }

    int endTime;
    tmr :> endTime;
    int timeTaken = endTime - startTime;
    printf("Time taken: %fms", timeTaken / 100000.0f);
}

void distributor(chanend c_in, chanend c_out, chanend fromAcc)
{
  uchar imageCurrent[IMHT][IMWD/8];
  uchar imageNext[IMHT][IMWD/8];

  //Starting up and wait for tilting of the xCore-200 Explorer
  printf( "ProcessImage: Start, size = %dx%d\n", IMHT, IMWD );
  printf( "Waiting for Board Tilt...\n" );
  fromAcc :> int value;

  printf( "Processing...\n" );
  for( int y = 0; y < IMHT; y++ ) {     //go through all lines
      for( int x = 0; x < IMWD / 8; x++ ) { //go through each pixel per line
        uchar nextByte = 0;
        for (int b = 0; b < 8; b++)
        {
            uchar tmp = 0;
            c_in :> tmp;
            if (tmp)
            {
                nextByte |= (1 << b);
            }
        }
        imageCurrent[y][x] = nextByte;
      }
  }

  generate(10, imageCurrent, imageNext);

  for( int y = 0; y < IMHT; y++ ) {   //go through all lines
    for( int x = 0; x < IMWD; x++ ) { //go through each pixel per line
      if (getCell(imageNext, x, y)) {
          c_out <: (uchar)255;
      }
      else {
          c_out <: (uchar)0;
      }
    }
  }
  printf( "\nOne processing round completed...\n" );
}

/////////////////////////////////////////////////////////////////////////////////////////
//
// Write pixel stream from channel c_in to PGM image file
//
/////////////////////////////////////////////////////////////////////////////////////////
void DataOutStream(char outfname[], chanend c_in)
{
  int res;
  uchar line[ IMWD ];

  //Open PGM file
  printf( "DataOutStream: Start...\n" );
  res = _openoutpgm( outfname, IMWD, IMHT );
  if( res ) {
    printf( "DataOutStream: Error opening %s\n.", outfname );
    return;
  }

  //Compile each line of the image and write the image line-by-line
  for( int y = 0; y < IMHT; y++ ) {
    for( int x = 0; x < IMWD; x++ ) {
      c_in :> line[ x ];
      //printf( "-%4.1d ", line[ x ] ); //show image values
    }
    _writeoutline( line, IMWD );
    printf( " DataOutStream: Line written...\n" );
  }

  //Close the PGM image
  _closeoutpgm();
  printf( "DataOutStream: Done...\n" );
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////
//
// Initialise and  read orientation, send first tilt event to channel
//
/////////////////////////////////////////////////////////////////////////////////////////
void orientation( client interface i2c_master_if i2c, chanend toDist) {
  i2c_regop_res_t result;
  char status_data = 0;
  int tilted = 0;

  // Configure FXOS8700EQ
  result = i2c.write_reg(FXOS8700EQ_I2C_ADDR, FXOS8700EQ_XYZ_DATA_CFG_REG, 0x01);
  if (result != I2C_REGOP_SUCCESS) {
    printf("I2C write reg failed\n");
  }
  
  // Enable FXOS8700EQ
  result = i2c.write_reg(FXOS8700EQ_I2C_ADDR, FXOS8700EQ_CTRL_REG_1, 0x01);
  if (result != I2C_REGOP_SUCCESS) {
    printf("I2C write reg failed\n");
  }

  //Probe the orientation x-axis forever
  while (1) {

    //check until new orientation data is available
    do {
      status_data = i2c.read_reg(FXOS8700EQ_I2C_ADDR, FXOS8700EQ_DR_STATUS, result);
    } while (!status_data & 0x08);

    //get new x-axis tilt value
    int x = read_acceleration(i2c, FXOS8700EQ_OUT_X_MSB);

    //send signal to distributor after first tilt
    if (!tilted) {
      if (x>30) {
        tilted = 1 - tilted;
        toDist <: 1;
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////
//
// Orchestrate concurrent system and start up all threads
//
/////////////////////////////////////////////////////////////////////////////////////////
int main(void) {

i2c_master_if i2c[1];               //interface to orientation

char infname[] = "test.pgm";     //put your input image path here
char outfname[] = "testout.pgm"; //put your output image path here
chan c_inIO, c_outIO, c_control;    //extend your channel definitions here

par {
    i2c_master(i2c, 1, p_scl, p_sda, 10);   //server thread providing orientation data
    orientation(i2c[0],c_control);        //client thread reading orientation data
    DataInStream(infname, c_inIO);          //thread to read in a PGM image
    DataOutStream(outfname, c_outIO);       //thread to write out a PGM image
    distributor(c_inIO, c_outIO, c_control);//thread to coordinate work on image
    //keepTime();
  }

  return 0;
}
