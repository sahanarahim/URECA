# Plant Connectome

This is a forked repository (i.e., repo) of the original plant connectome project.  This fork contains changes that I (i.e., Kevin) have made to the original repository.

Plant Connectome is a GPT3-powered database that is built on over 100,000 research paper abstracts and built with a Flask backend.  As of the time of writing, this project is still a work in progress with several tasks to be done.

More information about the project's files and folders can be found below (at least the major ones):

## Application Structure 

1.  `dics` (folder)

    This folder contains serialized Python objects (dictionaries mostly) that are crucial for Plant Connectome to run.  All data is currently stored in a seralized file `allDic` (or at least one of its variants).

1.  `infos` (folder)

    This folder contains text files that will be accessed during `app.py`'s runtime.

1.  `static` (folder)

    This folder contains images for now (though it can contain publicly-accessible attributes).

1.  `templates` (folder)

    Flask requires this folder to run for templating.  Each file inside this folder is a HTML file.

1.  `tests` (folder)

    This folder contains unit (and maybe integration tests in the future) tests for the application.

1.  `utils` (folder)

    This folder contains helper functions crucial for the application to run.

1.  `app.py` (file)

    This is the main application file.  Don't delete this one!

1.  `Procfile` (file)

    Heroku requires this file to run.  Don't delete this one as well!

1.  `requirements.txt` (file)

    This file contains a list of dependencies for this project.

## Setting up Plant Connectome

As Plant Connectome is a Flask application, Plant Connectome can also be set up in a typical manner as most Flask applications:

1.  First, ensure that you have Flask installed onto your machine.  If not, do install it via `pip` or similar software!
1.  Then, navigate to the project's root directory via the command line and set the Flask application via `set FLASK_APP=app.py`.  If you are using MacOS, type `export FLASK_APP=app.py` instead (i.e., replace `set` with `export`)!
1.  Then, type `flask run`.  A development server should begin, hence allowing you to see Plant Connectome in action on your machine!



