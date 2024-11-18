import numpy as np

# simulation parameters
SEED=0
VN_100 = 1
VN_002 = 1
ROWS=2
COLS=2
N_GRAINS=ROWS*COLS
SIZE=(100,100)
TIME_LIMIT=1000
STEP_TIME=1e-6
MAX_STEPS=1000
EPS=1e-10

# hexagonal grain constants
SQRT3 = np.sqrt(3)
TYPE_100 = 0
TYPE_002 = 1

ORIENTATION={   0:(0,1,0), 
                1:(SQRT3/2, 1/2, 0),
                2:(SQRT3/2, -1/2, 0),
                3:(0, -1, 0),
                4:(-SQRT3/2, -1/2, 0),
                5:(-SQRT3/2, 1/2, 0),
                6:(0,0,1),
                7:(0,0,-1)}

ADJACENCY={     0:[1,6,5,7],
                1:[2,6,0,7],
                2:[3,6,1,7],
                3:[4,6,2,7],
                4:[5,6,3,7],
                5:[0,6,4,7],
                6:[0,1,2,3,4,5],
                7:[0,1,2,3,4,5]}

TYPE={          0:TYPE_100,
                1:TYPE_100,
                2:TYPE_100,
                3:TYPE_100,
                4:TYPE_100,
                5:TYPE_100,
                6:TYPE_002,
                7:TYPE_002}

VN={            TYPE_100:VN_100,
                TYPE_002:VN_002}