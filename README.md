# BT-PPQ
A Burst-tolerant System to execute Parallel Preference Queries on Multicore Architectures

This framework supports the parallel execution of preference queries (e.g., skyline, top-k) in multicore systems using a sliding-window approach on data streams. The framework deals with bursty and time-varying input rates using a clever distribution strategies of input items to parallel entities within the system.

The code has been used in the paper "Parallel Continuous Preference Queries over Out-of-Order and Bursty Data Streams". IEEE Transactions on Parallel and Distributed Systems, 28(9), 2608-2624, 2017, IEEE.

The code requires FastFlow 2.1 (http://calvados.di.unipi.it/).
