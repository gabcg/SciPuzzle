#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Exceptions:
# Stoichiometry


class IncorrectStoichiometry(ValueError):
    """Exception for an Stoichiometry impossible to be done
    with the files given by the user."""

    def __init__(self, stoichiometry):
        self.stoichiometry = stoichiometry

    def __str__(self):
        a = 'The'
        b = 'can not be used.'
        return "%s %s %s\n" % (a, self.stoichiometry, b)
