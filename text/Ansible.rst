============================
Ansible: How I set up the VM
============================


Nobody expects the extra section

.. image:: imgs/monty-python-spanish-inquisition.png

-----------------------


`Ansible <https://docs.ansible.com/>`_ is an automation tool for setup of
virtual machines. It allows one to reliably and reproducibly perform setup
commands over a number of virtual machines.


It uses YAML files to define tasks. I'll live demo a setup of a vm that I use
for my section. It installs all required software and does a bit of
housekeeping that's worth doing.

I use ansible a lot for this sort of task. I can also show you how it can be
used to set up a docker image the same as a remove vm for debugging / writing
on the plane.
