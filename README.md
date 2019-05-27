# BT-PPQ
A Burst-Tolerant System to Execute Parallel Preference Queries on Multi-Core Architectures

This framework supports the parallel execution of preference queries (e.g., skyline, top-k) in multicore systems using a sliding-window approach on data streams. The framework deals with bursty and time-varying input rates using a clever distribution strategies of input items to parallel entities within the system.

The code requires FastFlow 2.1 (http://calvados.di.unipi.it/).

# Contributors
BT-PPQ has been developed by [Gabriele Mencagli](mailto:mencagli@di.unipi.it).

If our work is useful for your research, please cite the following paper:
 - G. Mencagli, M. Torquati, M. Danelutto and T. De Matteis. Parallel Continuous Preference Queries over Out-of-Order and Bursty Data Streams. IEEE Transactions on Parallel and Distributed Systems, 28(9), 2608-2624, 2017, IEEE. ISSN: 1045-9219, DOI: 10.1109/TPDS.2017.2679197
