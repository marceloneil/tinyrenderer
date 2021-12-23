#!/bin/sh

pqiv --watch-files=changes-only --action='set_interpolation_quality(0)' --action='toggle_scale_mode(3)' $1 &
