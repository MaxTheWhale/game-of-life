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
#define IMWD 64
#define IMHT 64
// Set the image you want to load in here, or generate one on-board
#define GENERATE 0
const char fname[] = "64x64.pgm";
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

// Loads the PGM file in and sends it down a channel
void DataInStream(const char infname[], chanend c_out)
{
    int res;
        uchar line[IMWD];
        printf("DataInStream: Start...\n");

        // If we're generating the image, then we make it a grid of 3-cell blinkers. This means the resulting image will have a lot of activity.
        if (GENERATE) {
            for (int y = 0; y < IMHT; y++) {
                for (int x = 0; x < IMWD; x++) {
                    // This if statement makes it so this pattern is repeated over and over
                    if (!(y % 4) && (x % 4)) {    // X X X .
                        c_out <: (uchar)255;      // . . . .
                    }                             // . . . .
                    else {                        // . . . .
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

// Gets the value of a cell given its x and y coordinates within a SLICE
uchar getCell(uchar slice[SLICEHT][LINE_SIZE], uint x, uint y){
    return slice[y][x / 8] & (1 << (x % 8)) ? 1 : 0;
}

// Processes the game of life for one line of a SLICE, writing the output into a separate line buffer
void processCells(uchar slice[SLICEHT][LINE_SIZE], uchar line[LINE_SIZE], uint y) {
    uint ym = (y - 1);
    uint yp = (y + 1);
    uint x = 0;
    uint xm = IMWD - 1;
    uint xp = 1;
    // Loops through every byte in the line
    for (int xb = 0; xb < LINE_SIZE; xb++) {
        uchar result = 0;
        // Calculates the new value for every bit in the byte
        for (int b = 0; b < 8; b++) {
            uchar living = getCell(slice, x, y);
            uchar neighbours = 0;

            // We count how many of the cells neighbours are alive
            neighbours += getCell(slice, xm, ym);
            neighbours += getCell(slice, x, ym);
            neighbours += getCell(slice, xp, ym);

            neighbours += getCell(slice, xm, y);
            neighbours += getCell(slice, xp, y);

            neighbours += getCell(slice, xm, yp);
            neighbours += getCell(slice, x, yp);
            neighbours += getCell(slice, xp, yp);

            // If the cell should be alive next round, we set the appropriate bit in the output
            if (neighbours == 3) {
                result |= (1 << b);
            }
            else if (neighbours == 2 && living) {
                result |= (1 << b);
            }
            // Increment coordinates for next cell
            x++;
            xm++; xm %= IMWD;
            xp++; xp %= IMWD;
        }
        // Write new byte to line buffer
        line[xb] = result;
    }
}

// Counts the number of currently living cells in a ROW
int liveCells(uchar slice[SLICEHT][LINE_SIZE]) {
    int count = 0;
    for (int y = 1; y < ROWHT + 1; y++) {
        for (int x = 0; x < IMWD; x++) {
            if (getCell(slice, x, y)) count++;
        }
    }
    return count;
}

// The main worker function
void worker(chanend c_dist) {
    uchar slice[SLICEHT][LINE_SIZE];
    uchar lineBuffer[2][LINE_SIZE];

    // First receive a ROW from the distributor
    for (int y = 1; y < ROWHT + 1; y++) {
        for (int x = 0; x < LINE_SIZE; x++) {
            c_dist :> slice[y][x];
        }
    }

    // Main loop begins
    while (1) {
        // Send the top and bottom lines of the ROW to the distributor
        for (int x = 0; x < LINE_SIZE; x++)
        {
            c_dist <: slice[1][x];
            c_dist <: slice[ROWHT][x];
        }
        // Receive the top and bottom lines of the SLICE back from the distributor
        for (int x = 0; x < LINE_SIZE; x++)
        {
            c_dist :> slice[0][x];
            c_dist :> slice[ROWHT + 1][x];
        }

        // Main loop to calculate the new ROW
        for (int y = 1; y < ROWHT + 1; y++) {
            // Write the new value of the line to the lineBuffer
            processCells(slice, lineBuffer[y % 2], y);
            for (int x = 0; x < LINE_SIZE; x++) {
                // Write back the previous line to the ROW
                if (y > 1) slice[y - 1][x] = lineBuffer[!(y % 2)][x];
                // If it is the final line then write it back immediately
                if (y == ROWHT) slice[y][x] = lineBuffer[y % 2][x];
            }
        }
        // Notify the distributor that we have finished processing
        c_dist <: (uchar)1;

        // Distributor tells us whether to export the image or not
        uchar exporting; c_dist :> exporting;
        if (exporting) {
            // Send our ROW back to the distributor
            for (int y = 1; y < ROWHT + 1; y++) {
                for (int x = 0; x < LINE_SIZE; x++) {
                    c_dist <: slice[y][x];
                }
            }
        }

        // Distributor tells us whether to send info or not
        uchar sendingInfo; c_dist :> sendingInfo;
        if (sendingInfo) {
            // Send the distributor the amount of live cells in our ROW
            c_dist <: liveCells(slice);
        }
        // Do not start the next generation until the distributor tells us to resume
        uchar resume; c_dist :> resume;
    }
}

// This function is used to receive the image from DataIn and send it to the workers
void distributeImage(chanend workers[WORKERS], chanend c_in) {
    for (int y = 0; y < IMHT; y++) {
        for (int x = 0; x < LINE_SIZE; x++) {
            // This loop packs the bytes we receive from DataIn into individual bits
            uchar nextByte = 0;
            for (int b = 0; b < 8; b++) {
                uchar tmp = 0;
                c_in :> tmp;
                if (tmp) {
                    nextByte |= (1 << b);
                }
            }
            // Send the packed byte to the appropriate worker
            workers[y / ROWHT] <: nextByte;
        }
    }
}

// Updates the tops and bottoms of the workers SLICEs ready for the next round
void updateSlices(chanend workers[WORKERS]) {
    // These are buffers to store the overlapping SLICE lines that the workers need from each other.
    uchar sliceTops[WORKERS][LINE_SIZE];
    uchar sliceBottoms[WORKERS][LINE_SIZE];
    // Receive the SLICE tops and bottoms from the workers
    for (int w = 0; w < WORKERS; w++) {
        for (int x = 0; x < LINE_SIZE; x++)
        {
            // The lines are sent to the appropriate place in the array with some arithmetic
            workers[w] :> sliceBottoms[(w + 7) % WORKERS][x];
            workers[w] :> sliceTops[(w + 1) % WORKERS][x];
        }
    }
    // Send the SLICE tops and bottoms to the workers
    for (int w = 0; w < WORKERS; w++) {
        for (int x = 0; x < LINE_SIZE; x++)
        {
            workers[w] <: sliceTops[w][x];
            workers[w] <: sliceBottoms[w][x];
        }
    }
}

// Exports the image by receiving it from the workers and sending it to DataOut
void exportImage(chanend workers[WORKERS], chanend c_out) {
    for (int y = 0; y < IMHT; y++) {
        for (int x = 0; x < LINE_SIZE; x++) {
            uchar nextByte = 0;
            // We receive the packed bytes from the workers and unpack them to send to DataOut
            workers[y / ROWHT] :> nextByte;
            for (int b = 0; b < 8; b++) {
                if (nextByte & (1 << b)) {
                    c_out <: (uchar)255;
                }
                else c_out <: (uchar)0;
            }
        }
    }
}

// Helper function to notify all workers
void notifyWorkers(uchar value, chanend workers[WORKERS]) {
    for (int w = 0; w < WORKERS; w++) {
        workers[w] <: value;
    }
}

// The distributor function is responsible for orchestrating the worker threads
void distributor(client interface time_iface time, client interface event_iface event, chanend c_in, chanend c_out, chanend workers[WORKERS]) {

    printf("ProcessImage: Start, size = %dx%d\n", IMWD, IMHT);
    printf("Waiting for SW1 press...\n");

    // Wait for the event handler to tell us to start
    select {
        case event.start(): break;
    }
    printf("Processing...\n");

    // Read in the image from DataInStream, send to workers
    distributeImage(workers, c_in);

    event.clear();
    delay_milliseconds(1000);
    int iterations = 0;
    int sTime = time.getTime();
    
    // Main loop begins
    while (1) {

        updateSlices(workers);
        // Wait for all workers to finish this round
        for (int w = 0; w < WORKERS; w++) {
            uchar waiting; workers[w] :> waiting;
        }
        iterations++;

        // Check if we've been told to export the image
        if (event.export()) {
            time.pause();
            // Notify the workers that we're exporting
            notifyWorkers(1, workers);

            // Send the current generation and the image to DataOutStream
            printf("\nExporting current image\n");
            c_out <: iterations;
            exportImage(workers, c_out);

            event.clear();
            time.unpause();
        }
        else {
            // Notify workers not to export
            notifyWorkers(0, workers);
        }

        // Check if we've been told to pause and print info
        if (event.sendInfo()) {
            time.pause();
            int totalLiving = 0;
            // Tell workers to send info and then receive their live cell counts
            notifyWorkers(1, workers);
            for (int w = 0; w < WORKERS; w++) {
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
            // Notify workers not to send info
            notifyWorkers(0, workers);
        }
        // Tell the event handler a new generation is about to begin
        event.newGeneration();
        // Notify workers to begin the next generation
        notifyWorkers(1, workers);
    }
}

// Writes out the current image. The file name has the current generation appended to it.
void DataOutStream(const char outfname[], chanend c_in) {
    int res;
    uchar line[IMWD];
    char file_name[40];

    // Loop forever so that we can export the image multiple times
    while (1) {
        // Receive the generation number from the distributor
        int currentGen;
        c_in :> currentGen;
        printf("DataOutStream: Start...\n");
        if (GENERATE) {
            // If we're generating the image, name the output based on the dimensions and current generation
            sprintf(file_name, "%dx%d_gen%d.pgm", IMWD, IMHT, currentGen);
        }
        else {
            // Otherwise, use the input file name and current generation
            sprintf(file_name, "%s_gen%d.pgm", outfname, currentGen);
        }

        // Open the PGM image for writing
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
    // Flags and LED state
    uchar sendInfo = 0;
    uchar exportImage = 0;
    int paused = 0;
    int started = 0;
    int leds = 0;

    // Timing variables
    int period = 10000000; // 100ms period
    timer tmr;
    unsigned time = 0;
    tmr :> time;

    while (1) {
        select {
            // We only poll the tilt and buttons every 100ms to be more efficient
            case tmr when timerafter(time + period) :> void:
                time += period;
                // Read from orientation thread and button port
                int tilt; fromTilt :> tilt;
                int buttons; buttonPort :> buttons;
                uchar sw1 = (buttons & 1) ? 0 : 1;
                uchar sw2 = (buttons & 2) ? 0 : 1;
                if (started) {
                    // If SW2 is pressed, set the export flag
                    if (sw2 && !exportImage) exportImage = 1;
                    // If we're not paused and tilted, set the send info flag
                    if (tilt && !sendInfo && !paused) sendInfo = 1;
                    // If we are paused and no longer tilted, unpause
                    if (!tilt && paused) {
                        i.start();
                    }
                }
                else if (sw1) {
                    // If SW1 is pressed, light the green LED and tell the distributor to start
                    leds |= 4;
                    ledPort <: leds;
                    started = 1;
                    paused = 1;
                    i.start();
                }
                break;
            case i.export() -> int output:
                // If we're exporting the image, light the blue LED and pause
                if (exportImage) {
                    leds |= 2;
                    ledPort <: leds;
                    paused = 1;
                    output = 1;
                }
                else output = 0;
                break;
            case i.sendInfo() -> int output:
                // If we're sending info, light the red LED and pause
                if (sendInfo) {
                    leds |= 8;
                    ledPort <: leds;
                    paused = 1;
                    output = 1;
                }
                else output = 0;
                break;
            case i.newGeneration():
                // Every generation, toggle the green LED and clear the flags
                leds ^= 1;
                ledPort <: leds;
                sendInfo = 0;
                exportImage = 0;
                break;
            case i.clear():
                // Once an event has finished, clear the RGB LED and pause flag
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

        // Run the first half the workers on tile 0
        par (int i = 0; i < WORKERS / 2; i++) {
            on tile[0]: worker(c_workers[i]);
        }
        // Run the second half of the workers on tile 1
        par (int i = WORKERS / 2; i < WORKERS; i++) {
            on tile[1]: worker(c_workers[i]);
        }
    }
    return 0;
}
