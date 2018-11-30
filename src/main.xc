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
#define IMHT 512
#define IMWD 512

// Various constants to help work with the image structure
#define LINE_SIZE (IMWD/8)
#define IMG_SIZE (IMHT*LINE_SIZE)
#define ROWHT (IMHT/WORKERS)
#define SLICEHT (ROWHT+2)
#define ROW_SIZE (ROWHT*LINE_SIZE)
#define SLICE_SIZE (SLICEHT*LINE_SIZE)

// The number of worker threads to run
#define WORKERS 8

typedef unsigned char uchar;

// Ports for the buttons, leds and i2c
on tile[0]: port buttons = XS1_PORT_4E;
on tile[0]: port leds = XS1_PORT_4F;
on tile[0]: port p_scl = XS1_PORT_1E;
on tile[0]: port p_sda = XS1_PORT_1F;

// Set the image you want to load in here
const char fname[] = "512x512.pgm";

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
    printf("DataInStream: Done...\n");
}

// Gets the value of a cell given its x and y coordinates
uchar getCell(uchar image[SLICEHT][LINE_SIZE], uchar x, uchar y){
    return image[y][x/8] & (1 << (x % 8)) ? 1 : 0;
}

// Sends a SLICE to worker i
void sendSlice(uchar image[IMHT][LINE_SIZE], int i, chanend c) {
    for (int y = (ROWHT * i) - 1; y < ((ROWHT * i) - 1 + SLICEHT); y++) {
        int iy = (y+IMHT) % IMHT;
        for (int x = 0; x < (LINE_SIZE); x++) {
            c <: image[iy][x];
        }
    }
}

// Set the cell at coordinates x, y based on the value v
void setCell(uchar image[ROWHT][LINE_SIZE], uchar x, uchar y, uchar v){
    if (v) {
        image[y][(x / 8)] |= (1 << (x % 8));
    }
    else {
        image[y][(x / 8)] &= ~(1 << (x % 8));
    }
}

// The main logic function for the game of life.
// It returns whether the cell at x, y should be alive or dead next round
// by counting the number of living neighbours.
uchar checkCell(uchar image[SLICEHT][LINE_SIZE], uchar x, uchar y) {
    uchar living = getCell(image, x, y);
    uchar neighbours = 0;
    uchar result = 0;
    uchar ym = (y - 1); ym %= IMHT;
    uchar yp = (y + 1); yp %= IMHT;
    uchar xm = (x - 1); xm %= IMWD;
    uchar xp = (x + 1); xp %= IMWD;

    neighbours += getCell(image, xm, ym);
    neighbours += getCell(image, x, ym);
    neighbours += getCell(image, xp, ym);

    neighbours += getCell(image, xm, y);
    neighbours += getCell(image, xp, y);

    neighbours += getCell(image, xm, yp);
    neighbours += getCell(image, x, yp);
    neighbours += getCell(image, xp, yp);

    if (living != 0) {
        result = (neighbours == 2 || neighbours == 3) ? 255 : 0;
    }
    else {
        result = (neighbours == 3) ? 255 : 0;
    }
    return result;
}

// The worker thread.
// It takes a SLICE from the distributor, and calculates the
// next ROW from it, which is then returned to the distributor.
void worker(chanend c) {
    uchar slice[SLICEHT][LINE_SIZE];
    uchar row[ROWHT][LINE_SIZE];
    while (1) {
        for (int y = 0; y < SLICEHT; y++) {
            for (int x = 0; x < LINE_SIZE; x++) {
                c :> slice[y][x];
            }
        }
        for (int y = 0; y < ROWHT; y++) {
            for (int x = 0; x < IMWD; x++) {
                setCell(row, x, y, checkCell(slice, x, (y + 1)));
            }
        }
        for (int y = 0; y < ROWHT; y++) {
            for (int x = 0; x < LINE_SIZE; x++) {
                c <: row[y][x];
            }
        }
    }
}

// Function to receive the ROW back from the worker
void getRow(uchar image[IMHT][LINE_SIZE], int i, chanend c) {
    for (int y = 0; y < (ROWHT); y++) {
        int iy = y + (ROWHT * i);
        for (int x = 0; x < (LINE_SIZE); x++) {
            c :> image[iy][x];
        }
    }
}

uchar getCellFinal(uchar image[IMHT][LINE_SIZE], uchar x, uchar y) {
    return image[y][x / 8] & (1 << (x % 8)) ? 1 : 0;
}

// Testing function that prints the current image to the console
void printWorld(uchar image[IMHT][LINE_SIZE]) {
    for (int y = 0; y < IMHT; y++) {
        for (int x = 0; x < IMWD; x++) {
            char p = (getCellFinal(image, x, y)) ? 'X' : '.';
            printf("%c", p);
        }
        printf("\n");
    }
}

// Counts the number of currently living cells
int liveCells(uchar image[IMHT][LINE_SIZE]) {
    int count = 0;
    for (int y = 0; y < IMHT; y++) {
        for (int x = 0; x < IMWD; x++) {
            if (getCellFinal(image, x, y)) count++;
        }
    }
    return count;
}

// The distributor function is responsible for farming out the sections of the image to the worker threads,
// and collecting the results back again. This process repeats continuously.
void distributor(client interface time_iface time, client interface event_iface event, chanend c_in, chanend c_out, chanend workers[WORKERS]) {
    uchar image[IMHT][LINE_SIZE];

    printf("ProcessImage: Start, size = %dx%d\n", IMWD, IMHT);
    printf("Waiting for SW1 press...\n");

    // Wait for the event handler to tell us to start
    select {
        case event.start(): break;
    }

    // Read in the image from DataInStream
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
            image[y][x] = nextByte;
        }
    }
    event.clear();

    int iterations = 0;
    int sTime = time.getTime();

    // Main loop begins
    while (1) {
        // Send SLICEs to each worker, and then receive ROWs back
        for (int i = 0; i < WORKERS; i++) {
            sendSlice(image, i, workers[i]);
        }
        for (int i = 0; i < WORKERS; i++) {
            getRow(image, i, workers[i]);
        }
        iterations++;

        // Check if we've been told to export the image
        if (event.export()) {
            time.pause();

            // Send the current generation and the image to DataOutStream
            printf("\nExporting current image\n");
            c_out <: iterations;
            for (int y = 0; y < IMHT; y++) {
                for (int x = 0; x < IMWD; x++) {
                    if (getCellFinal(image, x, y)) {
                        c_out <: (uchar)255;
                    }
                    else {
                        c_out <: (uchar)0;
                    }
                }
            }

            event.clear();
            time.unpause();
        }

        // Check if we've been told to pause and print info
        if (event.sendInfo()) {
            time.pause();

            printf("\n%d generations completed\n", iterations);
            printf("%d cells alive\n", liveCells(image));
            int totalTime = time.getTime() - sTime;
            printf("%dms total processing time, %.1fms per generation\n", totalTime, totalTime / (float)iterations);

            // Wait to be told to resume
            select {
                case event.start(): break;
            }
            event.clear();
            time.unpause();
        }

        // Inform the event handler that a new generation is beginning
        event.newGeneration();
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
        sprintf(file_name, "%s_gen%d.pgm", outfname, iteration);
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
        on tile[0]: {DataInStream(fname, c_inIO); DataOutStream(fname, c_outIO);} // Thread for image input and output (we can save a core as they don't need to run simultaneously)
        on tile[0]: distributor(time_if, event_if, c_inIO, c_outIO, c_workers);   // Thread to coordinate work on image
        par (int i = 0; i < WORKERS; i++) {
            on tile[1]: worker(c_workers[i]);                                     // Worker threads
        }
    }
    return 0;
}
