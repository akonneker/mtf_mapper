/*
  License: unknown
  Code based on jbkPava implementation by Jan de Leeuw
  Source: https://rpubs.com/deleeuw/262795
*/

#include "include/pava.h"
#include <vector>
using std::vector;

struct pava_block {
    double value;
    double weight;
    int size;
    int previous;
    int next;
};

void jbk_pava (double *x, int n) {
    vector<struct pava_block> blocks(n);
    
    for (int i = 0; i < n; i++) {
        blocks[i].value = x[i];
        blocks[i].weight = 1.0; // You could pass per-sample weights, but do we _have_ to?
        blocks[i].size = 1;
        blocks[i].previous = i - 1; 
        blocks[i].next = i + 1;   
    }    
    
    int active = 0;
    do {
        bool upsatisfied = false;
        int next = blocks[active].next; 
        if (next == n) {
            upsatisfied = true;
        } else {
            if (blocks[next].value > blocks[active].value) {
                upsatisfied = true;
            }
        }
          
        if (!upsatisfied) {
            double ww = blocks[active].weight + blocks[next].weight;
            int nextnext = blocks[next].next;
            double wxactive = blocks[active].weight * blocks[active].value;
            double wxnext = blocks[next].weight * blocks[next].value;
            blocks[active].value = (wxactive + wxnext) / ww;
            blocks[active].weight = ww;
            blocks[active].size += blocks[next].size;
            blocks[active].next = nextnext;
            if (nextnext < n) {
                blocks[nextnext].previous = active;
            }
            blocks[next].size = 0;
        }
        bool downsatisfied = false;
        int previous = blocks[active].previous;
        if (previous == -1) {
            downsatisfied = true;
        } else {
            if (blocks[previous].value < blocks[active].value) {
                downsatisfied = true;
            }
        }
        
        if (!downsatisfied) {
            double ww = blocks[active].weight + blocks[previous].weight;
            int previousprevious = blocks[previous].previous;
            double wxactive = blocks[active].weight * blocks[active].value;
            double wxprevious = blocks[previous].weight * blocks[previous].value;
            blocks[active].value = (wxactive + wxprevious) / ww;
            blocks[active].weight = ww;
            blocks[active].size += blocks[previous].size;
            blocks[active].previous = previousprevious;
            if (previousprevious > -1) {
                blocks[previousprevious].next = active;
            }
            blocks[previous].size = 0;
        }
        
        if ((blocks[active].next == n) && downsatisfied) {
            break;
        }
        
        if (upsatisfied && downsatisfied) {
            active = next;
        }
    } while (true);
    int k = 0;
    for (int i = 0; i < n; i++) {
        int blksize = blocks[i].size;
        if (blksize > 0.0) {
            for (int j = 0; j < blksize; j++) {
                x[k] = blocks[i].value;
                k++;
            }
        }
    }
}
  