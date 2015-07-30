#!/usr/bin/env python
# Author: Surhud More (surhud.more@ipmu.jp)
# v1.0 : 2015, July 29
import mechanize
import os
import getpass

'''
This python function uses the environmental variables HSC_USER and HSC_PASSWD to
authenticate against the HSC data server at NAOJ and submit a sql job. It should
be possible to modify the code to perform other tasks like generate cutouts,
etc.
'''
def submit_hsc_job(sql_job, br=None):
    br = mechanize.Browser()
    br.set_handle_robots(False)

    url = "https://hscdata.mtk.nao.ac.jp:4443/datasearch/"

    response = br.open(url)

    '''
    for i, form in enumerate(br.forms()):
        print "Form name:", i, form.name
        print form
    '''

    br.form = list(br.forms())[0]
    try:
        br['credential[account_name]'] = os.environ["HSC_USER"]
    except:
        br['credential[account_name]'] = raw_input(
            "Give your username for the HSC database: \n")

    try:
        br['credential[password]'] = os.environ["HSC_PASSWD"]
    except:
        br['credential[password]'] = getpass.getpass(
            "Give your password for the HSC database: \n")

    response = br.submit()

    '''
    for i, form in enumerate(br.forms()):
        print "Form name:", i, form.name
        print form
    '''


    br.form = list(br.forms())[0]
    br['catalog_job[sql]'] = sql_job

    response = br.submit(name='enqueue')

    # print response.read()
    print "Query submitted: Look for output here:"
    print "https://hscdata.mtk.nao.ac.jp:4443/datasearch/catalog_jobs"

    return br


if __name__ == "__main__":
    sql_job = "SELECT object_id  \
    FROM   ssp3_4_1_20141224.photoobj_mosaic__deepcoadd__iselect \
    WHERE  ra2000   BETWEEN  34.0 AND  36.0 \
    AND  decl2000 BETWEEN - 5.0 AND - 4.5 \
    AND  imag_kron < 25.5 \
    LIMIT 20"
    submit_hsc_job(sql_job)
