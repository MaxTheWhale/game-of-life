// COMS20001 - Cellular Automaton Farm
// (using the XMOS i2c accelerometer demo code)
//
// Throughout the code, we have used this terminology:
// A ROW represents the section of an image that an individual worker updates
// A SLICE represents the section of an image a worker requires to update the ROW
// (Meaning it has an extra line at the top and bottom)

#include <platform.h>
#include <xs1.h>
#include <stdio.h>
#include "pgmIO.h"
#include "i2c.h"

// Image width and height
// Width must be a multiple of 8
// Height must be a multiple of the number of workers
#define IMWD 512
#define IMHT 512
// Set the image you want to load in here, or generate one on-board
#define GENERATE 0
const char fname[] = "512x512.pgm";
// The number of worker threads to run
#define WORKERS 8

// Various constants to help work with the image structure
#define LINE_SIZE (IMWD/8)
#define IMG_SIZE (IMHT*LINE_SIZE)
#define ROWHT (IMHT/WORKERS)
#define SLICEHT (ROWHT+2)
#define ROW_SIZE (ROWHT*LINE_SIZE)
#define SLICE_SIZE (SLICEHT*LINE_SIZE)

typedef unsigned char uchar;

// Ports for the buttons, leds and i2c
on tile[0]: port buttons = XS1_PORT_4E;
on tile[0]: port leds = XS1_PORT_4F;
on tile[0]: port p_scl = XS1_PORT_1E;
on tile[0]: port p_sda = XS1_PORT_1F;

// Register addresses for orientation
#define FXOS8700EQ_I2C_ADDR 0x1E
#define FXOS8700EQ_XYZ_DATA_CFG_REG 0x0E
#define FXOS8700EQ_CTRL_REG_1 0x2A
#define FXOS8700EQ_DR_STATUS 0x0
#define FXOS8700EQ_OUT_X_MSB 0x1
#define FXOS8700EQ_OUT_X_LSB 0x2
#define FXOS8700EQ_OUT_Y_MSB 0x3
#define FXOS8700EQ_OUT_Y_LSB 0x4
#define FXOS8700EQ_OUT_Z_MSB 0x5
#define FXOS8700EQ_OUT_Z_LSB 0x6

// Interface for the timer
interface time_iface {
    int getTime();
    void pause();
    void unpause();
};

// Interface for the event handler
interface event_iface {
    [[notification]] slave void start();
    [[clears_notification]] void clear();
    int export();
    int sendInfo();
    void newGeneration();
};

// The timer function keeps track of the number of milliseconds elapsed
// It can be paused and resumed by the client
[[combinable]]
void keepTime(server interface time_iface i) {
    int period = 100000; // Gives 1ms time period
    int count = 0;
    int paused = 0;
    timer tmr;
    unsigned time = 0;
    tmr :> time;
    while (1) {
        select {
            case tmr when timerafter(time + period) :> void:
                if (!paused) count++;
                time += period;
                break;
            case i.getTime() -> int timeOut:
                timeOut = count;
                break;
            case i.pause():
                paused = 1;
                break;
            case i.unpause():
                paused = 0;
                break;
        }
    }
}

// Loads the pgm file in and sends it down a channel
void DataInStream(const char infname[], chanend c_out)
{
    int res;
        uchar line[IMWD];
        printf("DataInStream: Start...\n");

        if (GENERATE) {
            for (int y = 0; y < IMHT; y++) {
                for (int x = 0; x < IMWD; x++) {
                    if (!(y % 4) && (x % 4)) {
                        c_out <: (uchar)255;
                    }
                    else {
                        c_out <: (uchar)0;
                    }
                }
            }
        }
        else {
            // Open PGM file
            res = _openinpgm(infname, IMWD, IMHT);
            if (res) {
                printf("DataInStream: Error openening %s\n.", infname);
                return;
            }

            // Read image line-by-line and send byte by byte to channel c_out
            for (int y = 0; y < IMHT; y++) {
                _readinline(line, IMWD);
                for (int x = 0; x < IMWD; x++) {
                    c_out <: line[x];
                }
            }
            // Close PGM image file
            _closeinpgm();
        }
        printf("DataInStream: Done...\n");
}

// Gets the value of a cell given its x and y coordinates
uchar getCell(uchar slice[SLICEHT][LINE_SIZE], uint x, uint y){
    return slice[y][x / 8] & (1 << (x % 8)) ? 1 : 0;
}

// The main logic function for the game of life.
// It returns whether the cell at x, y should be alive or dead next round
// by counting the number of living neighbours.
uchar processCells(uchar slice[SLICEHT][LINE_SIZE], uint xb, uint y) {
    uchar result = 0;
    uint ym = (y - 1);
    uint yp = (y + 1);
    for (int b = 0; b < 8; b++) {
        uint x = (xb * 8) + b;
        uchar living = getCell(slice, x, y);
        uchar neighbours = 0;

        uint xm = (x - 1); xm %= IMWD;
        uint xp = (x + 1); xp %= IMWD;

        neighbours += getCell(slice, xm, ym);
        neighbours += getCell(slice, x, ym);
        neighbours += getCell(slice, xp, ym);

        neighbours += getCell(slice, xm, y);
        neighbours += getCell(slice, xp, y);

        neighbours += getCell(slice, xm, yp);
        neighbours += getCell(slice, x, yp);
        neighbours += getCell(slice, xp, yp);
        if (neighbours == 3) {
            result |= (1 << b);
        }
        else if (neighbours == 2 && living) {
            result |= (1 << b);
        }
    }
    return result;
}

// Counts the number of currently living cells
int liveCells(uchar slice[SLICEHT][LINE_SIZE]) {
    int count = 0;
    for (int y = 1; y < ROWHT + 1; y++) {
        for (int x = 0; x < IMWD; x++) {
            if (getCell(slice, x, y)) count++;
        }
    }
    return count;
}

// ODD workers receive top row from EVEN workers
// ODD workers receive bottom row from EVEN workers
// EVEN workers receive top row from ODD workers
// EVEN workers receive bottom row from ODD workers
void worker(chanend c_dist) {
    uchar slice[SLICEHT][LINE_SIZE];
    uchar lineBuffer[2][LINE_SIZE];
    for (int y = 1; y < ROWHT + 1; y++) {
        for (int x = 0; x < LINE_SIZE; x++) {
            c_dist :> slice[y][x];
        }
    }

    while (1) {
        for (int x = 0; x < LINE_SIZE; x++)
        {
            c_dist <: slice[1][x];
            c_dist <: slice[ROWHT][x];
        }
        for (int x = 0; x < LINE_SIZE; x++)
        {
            c_dist :> slice[0][x];
            c_dist :> slice[ROWHT + 1][x];
        }
        for (int y = 1; y < ROWHT + 1; y++) {
            for (int x = 0; x < LINE_SIZE; x++) {
                lineBuffer[y % 2][x] = processCells(slice, x, y);
            }
            for (int x = 0; x < LINE_SIZE; x++) {
                if (y > 1) slice[y - 1][x] = lineBuffer[!(y % 2)][x];
                if (y == ROWHT) slice[y][x] = lineBuffer[y % 2][x];
            }
        }
        c_dist <: (uchar)1;
        uchar exporting; c_dist :> exporting;
        if (exporting) {
            for (int y = 1; y < ROWHT + 1; y++) {
                for (int x = 0; x < LINE_SIZE; x++) {
                    c_dist <: slice[y][x];
                }
            }
        }
        uchar sendingInfo; c_dist :> sendingInfo;
        if (sendingInfo) {
            c_dist <: liveCells(slice);
        }
        uchar resume; c_dist :> resume;
    }
}

// The distributor function is responsible for farming out the sections of the image to the worker threads,
// and collecting the results back again. This process repeats continuously.
void distributor(client interface time_iface time, client interface event_iface event, chanend c_in, chanend c_out, chanend workers[WORKERS]) {

    printf("ProcessImage: Start, size = %dx%d\n", IMWD, IMHT);
    printf("Waiting for SW1 press...\n");

    // Wait for the event handler to tell us to start
    select {
        case event.start(): break;
    }

    // Read in the image from DataInStream, send to workers
    printf("Processing...\n");
    for (int y = 0; y < IMHT; y++) {
        for (int x = 0; x < LINE_SIZE; x++) {
            uchar nextByte = 0;
            for (int b = 0; b < 8; b++) {
                uchar tmp = 0;
                c_in :> tmp;
                if (tmp) {
                    nextByte |= (1 << b);
                }
            }
            workers[y / ROWHT] <: nextByte;
        }
    }
    event.clear();
    delay_milliseconds(5000);
    int iterations = 0;
    int sTime = time.getTime();
    uchar sliceTops[WORKERS][LINE_SIZE];
    uchar sliceBottoms[WORKERS][LINE_SIZE];
    
    // Main loop begins
    while (1) {
        for (int w = 0; w < WORKERS; w++) {
            for (int x = 0; x < LINE_SIZE; x++)
            {
                workers[w] :> sliceBottoms[(w + 7) % WORKERS][x];
                workers[w] :> sliceTops[(w + 1) % WORKERS][x];
            }
        }
        for (int w = 0; w < WORKERS; w++) {
            for (int x = 0; x < LINE_SIZE; x++)
            {
                workers[w] <: sliceTops[w][x];
                workers[w] <: sliceBottoms[w][x];
            }
        }
        for (int w = 0; w < WORKERS; w++) {
            uchar waiting; workers[w] :> waiting;
        }
        iterations++;

        // Check if we've been told to export the image
        if (event.export()) {
            time.pause();
            for (int w = 0; w < WORKERS; w++) {
                workers[w] <: (uchar)1;
            }

            // Send the current generation and the image to DataOutStream
            printf("\nExporting current image\n");
            c_out <: iterations;
            for (int y = 0; y < IMHT; y++) {
                for (int x = 0; x < LINE_SIZE; x++) {
                    uchar nextByte = 0;
                    workers[y / ROWHT] :> nextByte;
                    for (int b = 0; b < 8; b++) {
                        if (nextByte & (1 << b)) {
                            c_out <: (uchar)255;
                        }
                        else c_out <: (uchar)0;
                    }
                }
            }

            event.clear();
            time.unpause();
        }
        else {
            for (int w = 0; w < WORKERS; w++) {
                workers[w] <: (uchar)0;
            }
        }

        // Check if we've been told to pause and print info
        if (event.sendInfo() || iterations == 100) {
            time.pause();
            int totalLiving = 0;
            for (int w = 0; w < WORKERS; w++) {
                workers[w] <: (uchar)1;
                int living;
                workers[w] :> living;
                totalLiving += living;
            }

            printf("\n%d generations completed\n", iterations);
            printf("%d cells alive\n", totalLiving);
            int totalTime = time.getTime() - sTime;
            printf("%dms total processing time, %.1fms per generation\n", totalTime, totalTime / (float)iterations);

            // Wait to be told to resume
            select {
                case event.start(): break;
            }
            event.clear();
            time.unpause();
        }
        else {
            for (int w = 0; w < WORKERS; w++) {
                workers[w] <: (uchar)0;
            }
        }
        event.newGeneration();
        for (int w = 0; w < WORKERS; w++) {
            workers[w] <: (uchar)1;
        }
    }
}

// Writes out the current image. The file name has the current generation appended to it.
void DataOutStream(const char outfname[], chanend c_in) {
    int res;
    uchar line[IMWD];
    char file_name[40];

    while (1) {
        int iteration;
        c_in :> iteration;
        printf("DataOutStream: Start...\n");
        if (GENERATE) {
            sprintf(file_name, "%dx%d_gen%d.pgm", IMWD, IMHT, iteration);
        }
        else {
            sprintf(file_name, "%s_gen%d.pgm", outfname, iteration);
        }
        res = _openoutpgm(file_name, IMWD, IMHT);
        if (res) {
            printf("DataOutStream: Error opening %s\n.", file_name);
            return;
        }

        // Compile each line of the image and write the image line-by-line
        for (int y = 0; y < IMHT; y++) {
            for (int x = 0; x < IMWD; x++) {
                c_in :> line[x];
            }
            _writeoutline(line, IMWD);
            printf(" DataOutStream: Line written...\n");
        }

        // Close the PGM image
        _closeoutpgm();
        printf("DataOutStream: Done...\n");
    }
}

// Handles the reading of the orientation sensor via i2c
void orientation(client interface i2c_master_if i2c, chanend toEvent) {
    i2c_regop_res_t result;
    char status_data = 0;

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

    int tilted = 0;
    // Probe the orientation x-axis forever
    while (1) {

        // Check until new orientation data is available
        do {
            status_data = i2c.read_reg(FXOS8700EQ_I2C_ADDR, FXOS8700EQ_DR_STATUS, result);
        } while (!status_data & 0x08);

        // Get new x-axis tilt value
        int tilt = read_acceleration(i2c, FXOS8700EQ_OUT_X_MSB);
        if (tilt > 127) tilt = 0;

        // We use a different threshold once we are tilted so that it will not toggle sporadically
        if (tilted) {
            if (tilt > 25) {
                toEvent <: 1;
            }
            else {
                tilted = 0;
                toEvent <: 0;
            }
        }
        else {
            if (tilt > 30) {
                tilted = 1;
                toEvent <: 1;
            }
            else {
                toEvent <: 0;
            }
        }
    }
}

// This function handles all the LED and button I/O, as well as the tilt data.
// It controls whether the distributor should start, export the image, pause, or resume.
[[combinable]]
void eventHandler(in port buttonPort, out port ledPort, server interface event_iface i, chanend fromTilt) {
    uchar sendInfo = 0;
    uchar exportImage = 0;
    int paused = 0;
    int started = 0;
    int tilt = 0;
    int buttons = 0;
    int leds = 0;
    int period = 10000000; // 100ms period
    timer tmr;
    unsigned time = 0;
    tmr :> time;

    while (1) {
        select {
            // We only poll the tilt and buttons every 100ms to be more efficient
            case tmr when timerafter(time + period) :> void:
                time += period;
                fromTilt :> tilt;
                buttonPort :> buttons;
                uchar sw1 = (buttons & 1) ? 0 : 1;
                uchar sw2 = (buttons & 2) ? 0 : 1;
                if (started) {
                    if (sw2 && !exportImage) exportImage = 1;
                    if (tilt && !sendInfo && !paused) sendInfo = 1;
                    if (!tilt && paused) {
                        i.start();
                    }
                }
                else if (sw1) {
                    leds |= 4;
                    ledPort <: leds;
                    started = 1;
                    paused = 1;
                    i.start();
                }
                break;
            case i.export() -> int output:
                if (exportImage) {
                    leds |= 2;
                    ledPort <: leds;
                    paused = 1;
                    output = 1;
                }
                else output = 0;
                break;
            case i.sendInfo() -> int output:
                if (sendInfo) {
                    leds |= 8;
                    ledPort <: leds;
                    paused = 1;
                    output = 1;
                }
                else output = 0;
                break;
            case i.newGeneration():
                leds ^= 1;
                ledPort <: leds;
                sendInfo = 0;
                exportImage = 0;
                break;
            case i.clear():
                leds &= ~14;
                ledPort <: leds;
                paused = 0;
                break;
        }
    }
}

// In main we define our channels, and then set up all the threads of our concurrent system.
int main() {
    // Interface definitions
    i2c_master_if i2c[1];
    interface time_iface time_if;
    interface event_iface event_if;

    // Channel definitions
    chan c_inIO, c_outIO, c_control;
    chan c_workers[WORKERS];

    par {
        on tile[0].core[0]: keepTime(time_if);                                    // Timer thread
        on tile[0].core[0]: eventHandler(buttons, leds, event_if, c_control);     // Event handler
        on tile[0]: i2c_master(i2c, 1, p_scl, p_sda, 10);                         // Server thread providing orientation data
        on tile[0]: orientation(i2c[0], c_control);                               // Client thread reading orientation data
        on tile[1]: {DataInStream(fname, c_inIO); DataOutStream(fname, c_outIO);} // Thread for image input and output (we can save a core as they don't need to run simultaneously)
        on tile[1]: distributor(time_if, event_if, c_inIO, c_outIO, c_workers);   // Thread to coordinate work on image
        par (int i = 0; i < WORKERS / 2; i++) {
            on tile[0]: worker(c_workers[i]);                                     // Worker threads
        }
        par (int i = WORKERS / 2; i < WORKERS; i++) {
            on tile[1]: worker(c_workers[i]);                                     // Worker threads
        }
    }
    return 0;
}
