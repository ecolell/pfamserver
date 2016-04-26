# -*- coding: utf-8 -*-
from __future__ import print_function
import schedule
from threading import Thread, Event
import time
import os
import traceback


def run_schedule(stop_t):
    import core
    while not stop_t.isSet():
        schedule.run_pending()
        stop_t.wait(10)


def run_every(moment="day", hour="00:00"):

    def real_decorator(function):
        print(" * --> Schedule {:} to execute every {:} at {:}".
              format(function.__name__, moment, hour))
        getattr(schedule.every(), moment).at(hour).do(function)

        def wrapper(*args, **kwargs):
            try:
                p = Process(target=function, args=args, kwargs=kwargs)
                p.start()
            except:
                print(traceback.format_exc())
        return wrapper
    return real_decorator


def init_app(app):
    app.run_every = run_every
    if not app.debug or os.environ.get("WERKZEUG_RUN_MAIN") == "true":
        app.stop_t = Event()
        app.t = Thread(target=run_schedule, args=(app.stop_t, ))
        app.t.start()


def finish_app(app):
    if hasattr(app, 'stop_t'):
        app.stop_t.set()
    if hasattr(app, 't'):
        app.t.join()

