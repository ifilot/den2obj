.. _faq:
.. index:: Frequently Asked Questions

Frequently Asked Questions
==========================

.. contents::
   :local:

..
  Frequently asked questions should be questions that actually got asked.
  Formulate them as a question and an answer.
  Consider that the answer is best as a reference to another place in the documentation.


Creating .d2o files
-------------------

Why do I see the error "Decompression could not be verified."?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After compression, :program:`Den2Obj` checks whether the compression
is succesfull by decompressing the stream and tesing for similarity.
If the decompressed stream is different from the original stream, this
error is thrown automatically. Although it is very rare for this error
to be thrown, we have seen it happening for LZMA compression. The
recommended work-around is to use a different compression routine,
such as ``gzip`` or ``bzip2``.

.. seealso::

   :doc:`user_interface`
      Common errors and solutions for build failures.