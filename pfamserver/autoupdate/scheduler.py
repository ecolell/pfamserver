# -*- coding: utf-8 -*-
from __future__ import print_function
import schedule
from threading import Thread
from multiprocessing import Process
import time
import os
import traceback


def run_schedule():
    import core
    while 1:
        schedule.run_pending()
        time.sleep(10)


def run_every(moment="day", hour="00:00"):

    def real_decorator(function):
        print(" * --> Schedule {:} to execute every {:} at {:}".
              format(function.__name__, moment, hour))
        getattr(schedule.every(), moment).at(hour).do(function)

        def wrapper(*args, **kwargs):
            try:
                function(*args, **kwargs)
            except:
                print(traceback.format_exc())
        return wrapper
    return real_decorator


def init_app(app):
    app.run_every = run_every
    if not app.debug or os.environ.get("WERKZEUG_RUN_MAIN") == "true":
        t = Thread(target=run_schedule)
        t.start()
