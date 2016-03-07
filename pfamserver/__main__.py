from application import app
import os
from autoupdate import update
from loader import init_db
import os


if __name__ == '__main__':
    if not os.environ.get("WERKZEUG_RUN_MAIN"):
        #TODO: It download in a parallel process. And when reboot it
        #should install.
        pass
        # init_db()
        # update()
        print "Loading"
    import api
    import admin
    import web
    host = os.getenv('HOST', '0.0.0.0')
    port = os.getenv('PORT', 5001)
    app.run(host=host, port=port, threaded=True)
