from ..draw import bond_source

source = bond_source(8)
bond_source = {
    4: (
        [
            [0.0, 1.0, -0.5],
            [0.0, 1.0, 0.5],
            [1.0, -4.371138828673793e-08, -0.5],
            [1.0, -4.371138828673793e-08, 0.5],
            [-8.742277657347586e-08, -1.0, -0.5],
            [-8.742277657347586e-08, -1.0, 0.5],
            [-1.0, 1.1924880638503055e-08, -0.5],
            [-1.0, 1.1924880638503055e-08, 0.5],
        ],
        [
            [0, 1, 3, 2],
            [2, 3, 5, 4],
            [3, 1, 7, 5],
            [4, 5, 7, 6],
            [6, 7, 1, 0],
            [0, 2, 4, 6],
        ],
        [],
    ),
    6: (
        [
            [0.0, 1.0, -0.5],
            [0.0, 1.0, 0.5],
            [0.866025447845459, 0.4999999701976776, -0.5],
            [0.866025447845459, 0.4999999701976776, 0.5],
            [0.8660253882408142, -0.5000000596046448, -0.5],
            [0.8660253882408142, -0.5000000596046448, 0.5],
            [-8.742277657347586e-08, -1.0, -0.5],
            [-8.742277657347586e-08, -1.0, 0.5],
            [-0.866025447845459, -0.49999991059303284, -0.5],
            [-0.866025447845459, -0.49999991059303284, 0.5],
            [-0.866025447845459, 0.49999991059303284, -0.5],
            [-0.866025447845459, 0.49999991059303284, 0.5],
        ],
        [
            [0, 1, 3, 2],
            [2, 3, 5, 4],
            [4, 5, 7, 6],
            [6, 7, 9, 8],
            [8, 9, 11, 10],
            [10, 11, 1, 0],
        ],
        [[3, 1, 11, 9, 7, 5], [0, 2, 4, 6, 8, 10]],
    ),
    8: (
        [
            [0.0, 1.0, -0.5],
            [0.0, 1.0, 0.5],
            [0.7071067690849304, 0.7071067690849304, -0.5],
            [0.7071067690849304, 0.7071067690849304, 0.5],
            [1.0, -4.371138828673793e-08, -0.5],
            [1.0, -4.371138828673793e-08, 0.5],
            [0.7071067690849304, -0.7071067690849304, -0.5],
            [0.7071067690849304, -0.7071067690849304, 0.5],
            [-8.742277657347586e-08, -1.0, -0.5],
            [-8.742277657347586e-08, -1.0, 0.5],
            [-0.7071067094802856, -0.7071068286895752, -0.5],
            [-0.7071067094802856, -0.7071068286895752, 0.5],
            [-1.0, 1.1924880638503055e-08, -0.5],
            [-1.0, 1.1924880638503055e-08, 0.5],
            [-0.70710688829422, 0.7071066498756409, -0.5],
            [-0.70710688829422, 0.7071066498756409, 0.5],
        ],
        [
            [0, 1, 3, 2],
            [2, 3, 5, 4],
            [4, 5, 7, 6],
            [6, 7, 9, 8],
            [8, 9, 11, 10],
            [10, 11, 13, 12],
            [12, 13, 15, 14],
            [14, 15, 1, 0],
        ],
        [[3, 1, 15, 13, 11, 9, 7, 5], [0, 2, 4, 6, 8, 10, 12, 14]],
    ),
    12: (
        [
            [0.0, 1.0, -0.5],
            [0.0, 1.0, 0.5],
            [0.5, 0.8660253882408142, -0.5],
            [0.5, 0.8660253882408142, 0.5],
            [0.866025447845459, 0.4999999701976776, -0.5],
            [0.866025447845459, 0.4999999701976776, 0.5],
            [1.0, -4.371138828673793e-08, -0.5],
            [1.0, -4.371138828673793e-08, 0.5],
            [0.8660253882408142, -0.5000000596046448, -0.5],
            [0.8660253882408142, -0.5000000596046448, 0.5],
            [0.5000000596046448, -0.8660253882408142, -0.5],
            [0.5000000596046448, -0.8660253882408142, 0.5],
            [-8.742277657347586e-08, -1.0, -0.5],
            [-8.742277657347586e-08, -1.0, 0.5],
            [-0.4999999701976776, -0.8660253882408142, -0.5],
            [-0.4999999701976776, -0.8660253882408142, 0.5],
            [-0.866025447845459, -0.49999991059303284, -0.5],
            [-0.866025447845459, -0.49999991059303284, 0.5],
            [-1.0, 1.1924880638503055e-08, -0.5],
            [-1.0, 1.1924880638503055e-08, 0.5],
            [-0.866025447845459, 0.49999991059303284, -0.5],
            [-0.866025447845459, 0.49999991059303284, 0.5],
            [-0.5000001788139343, 0.8660253286361694, -0.5],
            [-0.5000001788139343, 0.8660253286361694, 0.5],
        ],
        [
            [0, 1, 3, 2],
            [2, 3, 5, 4],
            [4, 5, 7, 6],
            [6, 7, 9, 8],
            [8, 9, 11, 10],
            [10, 11, 13, 12],
            [12, 13, 15, 14],
            [14, 15, 17, 16],
            [16, 17, 19, 18],
            [18, 19, 21, 20],
            [20, 21, 23, 22],
            [22, 23, 1, 0],
        ],
        [
            [3, 1, 23, 21, 19, 17, 15, 13, 11, 9, 7, 5],
            [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22],
        ],
    ),
    16: (
        [
            [0.0, 1.0, -0.5],
            [0.0, 1.0, 0.5],
            [0.3826834559440613, 0.9238795042037964, -0.5],
            [0.3826834559440613, 0.9238795042037964, 0.5],
            [0.7071067690849304, 0.7071067690849304, -0.5],
            [0.7071067690849304, 0.7071067690849304, 0.5],
            [0.9238795042037964, 0.3826834261417389, -0.5],
            [0.9238795042037964, 0.3826834261417389, 0.5],
            [1.0, -4.371138828673793e-08, -0.5],
            [1.0, -4.371138828673793e-08, 0.5],
            [0.9238795638084412, -0.3826833963394165, -0.5],
            [0.9238795638084412, -0.3826833963394165, 0.5],
            [0.7071067690849304, -0.7071067690849304, -0.5],
            [0.7071067690849304, -0.7071067690849304, 0.5],
            [0.38268348574638367, -0.9238795042037964, -0.5],
            [0.38268348574638367, -0.9238795042037964, 0.5],
            [-8.742277657347586e-08, -1.0, -0.5],
            [-8.742277657347586e-08, -1.0, 0.5],
            [-0.3826834261417389, -0.9238795042037964, -0.5],
            [-0.3826834261417389, -0.9238795042037964, 0.5],
            [-0.7071067094802856, -0.7071068286895752, -0.5],
            [-0.7071067094802856, -0.7071068286895752, 0.5],
            [-0.9238795042037964, -0.38268357515335083, -0.5],
            [-0.9238795042037964, -0.38268357515335083, 0.5],
            [-1.0, 1.1924880638503055e-08, -0.5],
            [-1.0, 1.1924880638503055e-08, 0.5],
            [-0.9238794445991516, 0.3826836049556732, -0.5],
            [-0.9238794445991516, 0.3826836049556732, 0.5],
            [-0.70710688829422, 0.7071066498756409, -0.5],
            [-0.70710688829422, 0.7071066498756409, 0.5],
            [-0.3826834261417389, 0.9238795638084412, -0.5],
            [-0.3826834261417389, 0.9238795638084412, 0.5],
        ],
        [
            [0, 1, 3, 2],
            [2, 3, 5, 4],
            [4, 5, 7, 6],
            [6, 7, 9, 8],
            [8, 9, 11, 10],
            [10, 11, 13, 12],
            [12, 13, 15, 14],
            [14, 15, 17, 16],
            [16, 17, 19, 18],
            [18, 19, 21, 20],
            [20, 21, 23, 22],
            [22, 23, 25, 24],
            [24, 25, 27, 26],
            [26, 27, 29, 28],
            [28, 29, 31, 30],
            [30, 31, 1, 0],
        ],
        [
            [3, 1, 31, 29, 27, 25, 23, 21, 19, 17, 15, 13, 11, 9, 7, 5],
            [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30],
        ],
    ),
}
