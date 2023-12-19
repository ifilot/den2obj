.. _d2ofileformat:
.. index:: D2O file format

D2O File Format
===============

The ``.d2o`` file format is a binary file type to store densely packed scalar
fields. :program:`Den2Obj` uses this file format to store the data in the
relatively unwieldy ``CHGCAR`` and ``PARCHG``. In general, ``.d2o`` files are
about a factor 5 smaller than their corresponding ``CHGCAR`` files and can be
read in only a fraction of the time because the data does not need to be
converted from a human-readable format to binary.

The organization of this file is given in the Table 1.

.. list-table:: Table 1: Storage format of the ``.d2o`` file.
    :header-rows: 1
    :class: tight-table

    * - Position
      - Type
      - Length
      - Description
    * - 0x00-0x02
      - char
      - 3 bytes
      - Fixed token of "D2O" to identify the file.
    * - 0x03-0x07
      - uint32_t
      - 4 bytes
      - Protocol identifier token (see below). Default = 2.
    * - 0x08-0x2B
      - float[9]
      - 36 bytes
      - Unit cell matrix
    * - 0x2C-0x37
      - uint32_t[3]
      - 12 bytes
      - Grid dimensions (its product is the number of data points)
    * - 0x38
      - uint8_t
      - 1 byte
      - Floating point bytesize (float = 4, double = 8). Mainly used for validation purposes.
    * - 0x39-0x41
      - uint64_t
      - 8 bytes
      - Size of the compressed data stream
    * - 0x42..
      - char[DATASIZE]
      - DATASIZE bytes
      - Compressed data stream containing densely packed scalar field.

To convert a ``CHGCAR`` file to ``.d2o`` file format, run the following command::

    den2obj -i CHGCAR -o filename.d2o -t

.. note::
   :program:`Den2Obj` will automatically look for the best compression algorithm
   when converting a scalar field to the ``.d2o`` format. For the majority of the
   cases, this corresponds to the `LZMA type of compression <https://en.wikipedia.org/wiki/Lempel%E2%80%93Ziv%E2%80%93Markov_chain_algorithm>`_.

Protocol tokens
---------------

1. GZIP compression
2. LZMA compression
3. BZIP2 compression
