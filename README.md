# BT-PPQ
A Burst-Tolerant System to Execute Parallel Preference Queries on Multi-Core Architectures

This framework supports the parallel execution of preference queries (e.g., skyline, top-k) in multicore systems using a sliding-window approach on data streams. The framework deals with bursty and time-varying input rates using a clever distribution strategies of input items to parallel entities within the system.

The code requires FastFlow 2.1 (http://calvados.di.unipi.it/).

# How to Cite
If our work is useful for your research, please cite the following paper:
```
@ARTICLE{7873332,
  author={Mencagli, Gabriele and Torquati, Massimo and Danelutto, Marco and De Matteis, Tiziano},
  journal={IEEE Transactions on Parallel and Distributed Systems},
  title={Parallel Continuous Preference Queries over Out-of-Order and Bursty Data Streams},
  year={2017},
  volume={28},
  number={9},
  pages={2608-2624},
  doi={10.1109/TPDS.2017.2679197}
}
```

# Contributors
BT-PPQ has been developed by [Gabriele Mencagli](mailto:gabriele.mencagli@unipi.it)
