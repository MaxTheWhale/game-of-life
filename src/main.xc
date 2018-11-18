// COMS20001 - Cellular Automaton Farm - Initial Code Skeleton
// (using the XMOS i2c accelerometer demo code)

#include <platform.h>
#include <xs1.h>
#include <stdio.h>
#include "pgmIO.h"
#include "i2c.h"

#define IMHT 16                //image height
#define IMWD 16                //image width
#define LINE_SIZE (IMWD/8)
#define IMG_SIZE (IMHT*LINE_SIZE)
#define WORKERS 8
#define ROWHT (IMHT/WORKERS)
#define SLICEHT (ROWHT+2)
#define ROW_SIZE (ROWHT*LINE_SIZE)
#define SLICE_SIZE (SLICEHT*LINE_SIZE)

typedef unsigned char uchar;      //using uchar as shorthand
typedef unsigned long ulong;

on tile[0]: port p_scl = XS1_PORT_1E;         //interface ports to orientation
on tile[0]: port p_sda = XS1_PORT_1F;

char infname[] = "test.pgm";     //put your input image path here
char outfname[] = "testout.pgm"; //put your output image path here

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

interface time_iface {
    int getTime();
};

void keepTime(server interface time_iface i) {
  int count = 0;
  int period = 100000;
  timer tmr;
  unsigned time;
  while (1) {
    select {
      case tmr when timerafter(time) :> void:
        count++;
        time += period;
        break;
      case i.getTime() -> int timeOut:
        timeOut = count;
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

uchar getCell(uchar image[SLICEHT][LINE_SIZE], uchar x, uchar y){
    return image[y][x/8]&(1<<(x%8)) ? 1 : 0;
}

void printSlice(uchar slice[SLICEHT][LINE_SIZE]) {
    for( int y = 0; y < SLICEHT; y++ ) {
        for( int x = 0; x < IMWD; x++ ) {
          char p = (getCell(slice, x, y)) ? 'X' : '.';
          printf( "%c ",p); //show image values
        }
        printf( "\n" );
      }
}

void sendSlice(uchar image[IMHT][LINE_SIZE], int i, chanend c) {
    for (int y = (ROWHT * i) - 1; y < ((ROWHT * i) - 1 + SLICEHT); y++) {
        int iy = (y+IMHT) % IMHT;
        for (int x = 0; x < (LINE_SIZE); x++) {
          c <: image[iy][x];
        }
    }
}

void setCell(uchar image[ROWHT][LINE_SIZE], uchar x, uchar y, uchar v){
    if (v){
        image[y][(x/8)] |= (1<<(x%8));
    } else {
        image[y][(x/8)] &= ~(((uchar)1)<<(x%8));
    }
}

uchar checkCell(uchar image[SLICEHT][LINE_SIZE], uchar x, uchar y)
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

void worker(chanend c) {
  uchar slice[SLICEHT][LINE_SIZE];
  uchar row[ROWHT][LINE_SIZE];
  while(1) {
    for (int y = 0; y < SLICEHT; y++) {
      for (int x = 0; x < LINE_SIZE; x++) {
        c :> slice[y][x];
      }
    }
    /*printf("\n");
    printSlice(slice);*/
    //printf("Workers got here before the other statement\n");
    for (int y = 0; y < ROWHT; y++) {   //go through all lines
      for( int x = 0; x < IMWD; x++ ) {
        setCell(row, x, y, checkCell(slice, x, (y+1)));
      }
    }
    //printf("Workers got here\n");
    for (int y = 0; y < ROWHT; y++) {
      for (int x = 0; x < LINE_SIZE; x++) {
        c <: row[y][x];
      }
    }/*
    for (int y = 1; y < (IMHT/WORKERS+1); y++) {
      for (int x = 0; x < (IMWD/8); x++) {
        c <: slice[y][x];
      }
    }*/
  }
}

uchar getCellRow(uchar image[ROWHT][LINE_SIZE], uchar x, uchar y){
    return image[y][x/8]&(1<<(x%8)) ? 1 : 0;
}

void printRow(uchar image[ROWHT][LINE_SIZE]) {
    for( int y = 0; y < ROWHT; y++ ) {
        for( int x = 0; x < IMWD; x++ ) {
          char p = (getCellRow(image, x, y)) ? 'X' : '.';
          printf( "%c ",p); //show image values
        }
        printf( "\n" );
      }
}

void getRow(uchar image[IMHT][LINE_SIZE], int i, chanend c){
  for (int y = 0; y < (ROWHT); y++) {
    int iy = y + (ROWHT * i);
    for (int x = 0; x < (LINE_SIZE); x++) {
      c :> image[iy][x];
    }
  }
  /*int k = 0;
  for(int y = (IMHT/WORKERS) * i; y < ((IMHT/WORKERS) * (i+1)); y++){
      memcpy(image[y], r[k], sizeof(image[y]));
      k++;
  }*/
}

uchar getCellFinal(uchar image[IMHT][LINE_SIZE], uchar x, uchar y){
    return image[y][x/8]&(1<<(x%8)) ? 1 : 0;
}

void printWorld(uchar image[IMHT][LINE_SIZE]) {
    for( int y = 0; y < IMHT; y++ ) {
        for( int x = 0; x < IMWD; x++ ) {
          char p = (getCellFinal(image, x, y)) ? 'X' : '.';
          printf( "%c ",p); //show image values
        }
        printf( "\n" );
      }
}

void distributor(client interface time_iface i, chanend c_in, chanend c_out, chanend fromAcc, chanend workers[WORKERS])
{
  uchar image[IMHT][LINE_SIZE];

  //Starting up and wait for tilting of the xCore-200 Explorer
  printf( "ProcessImage: Start, size = %dx%d\n", IMHT, IMWD );
  printf( "Waiting for Board Tilt...\n" );
  fromAcc :> int value;

  printf( "Processing...\n" );
  for( int y = 0; y < IMHT; y++ ) {     //go through all lines
      for( int x = 0; x < LINE_SIZE; x++ ) { //go through each pixel per line
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
        image[y][x] = nextByte;
      }
  }
  int sTime = i.getTime();
  timer tmr;
  unsigned int startTime;
  tmr :> startTime;
  //printWorld(image);
  //printf("\n");
  for (int iteration = 0; iteration < 10; iteration++) {
      for (int i = 0; i < WORKERS; i++) {
          sendSlice(image, i, workers[i]);
        }
        //printf("got here\n");

        for (int i = 0; i < WORKERS; i++) {
          //printf("getting row %d\n",i);
          getRow(image, i, workers[i]);
        }
        //printf("got here, too\n");

        //printWorld(image);
        //printf("\n");
  }
  unsigned int endTime;
  tmr :> endTime;
  int timeTaken = endTime - startTime;
  int eTime = i.getTime();
  printf("Time taken: %.2fms  ", timeTaken / 100000.0f);
  printf("Time taken: %dms  ", (eTime - sTime));
  for( int y = 0; y < IMHT; y++ ) {   //go through all lines
    for( int x = 0; x < IMWD; x++ ) { //go through each pixel per line
      if (getCellFinal(image, x, y)) {
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
interface time_iface i;

chan c_inIO, c_outIO, c_control;  //extend your channel definitions here
chan workers[WORKERS];

par {
    on tile[0]: i2c_master(i2c, 1, p_scl, p_sda, 10);   //server thread providing orientation data
    on tile[0]: orientation(i2c[0],c_control);        //client thread reading orientation data
    on tile[0]: DataInStream(infname, c_inIO);          //thread to read in a PGM image
    on tile[0]: DataOutStream(outfname, c_outIO);       //thread to write out a PGM image
    on tile[0]: distributor(i, c_inIO, c_outIO, c_control, workers); //thread to coordinate work on image
    on tile[0]: keepTime(i);
    par (int i = 0; i < WORKERS; i++) {
        on tile[1]: worker(workers[i]);
    }
  }

  return 0;
}
