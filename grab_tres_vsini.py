import numpy as np
import pytest
import time
import json
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.support import expected_conditions
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.desired_capabilities import DesiredCapabilities
from io import StringIO

# Generated by Selenium IDE
class AccessingTresRVs():
  def setup_method(self, method):
    self.driver = webdriver.Chrome()
    self.vars = {}
  
  def teardown_method(self, method):
    self.driver.quit()
  
  def wait_for_window(self, timeout = 2):
    time.sleep(round(timeout / 1000))
    wh_now = self.driver.window_handles
    wh_then = self.vars["window_handles"]
    if len(wh_now) > len(wh_then):
      return set(wh_now).difference(set(wh_then)).pop()
  
  def get_table(self, usr, pwd, TICID):
    '''
    Scrapes TRES/CHIRON site to grab table
    '''

    TICID = str(TICID)
    TICID = TICID.zfill(9)
    T0name = 'T0' + TICID
    self.driver.get("http://tess.exoplanets.dk/Login.aspx")
    #self.driver.set_window_size(1397, 861)
    self.driver.find_element(By.ID, "ctl00_cpMainContent_tbUserName").click()
    self.driver.find_element(By.ID, "ctl00_cpMainContent_tbUserName").send_keys(usr)
    self.driver.find_element(By.ID, "ctl00_cpMainContent_tbPassword").send_keys(pwd)
    self.driver.find_element(By.ID, "ctl00_cpMainContent_btnLogin").click()
    self.driver.find_element(By.ID, "ctl00_cpMainContent_tbSearch").click()
    self.driver.find_element(By.ID, "ctl00_cpMainContent_tbSearch").send_keys(T0name)
    self.driver.find_element(By.ID, "ctl00_cpMainContent_tbSearch").send_keys(Keys.ENTER)
    self.driver.find_element(By.LINK_TEXT, T0name).click()
    self.vars["window_handles"] = self.driver.window_handles
    self.driver.find_element(By.ID, "ctl00_cpMainContent_Button12").click()
    self.vars["win6497"] = self.wait_for_window(2000)
    self.driver.switch_to.window(self.vars["win6497"])
    self.vars["fehtable"] = self.driver.find_element(By.CSS_SELECTOR, "pre").text
  
def grab_tres_vsini(username, password, TICID):
    '''
    Grabs vsini measurements and calculates the mean and standard error of the mean.
    '''

    myClass = AccessingTresRVs()
    myClass.setup_method("")
    myClass.get_table(usr=username, pwd=password, TICID=TICID)
    myClass.teardown_method("")

    datatable = StringIO(myClass.vars["fehtable"])

    # f = open(path + rvfilename, "w")
    # f.write(myClass.vars["rvtable"])
    # f.close()
    data = pd.read_csv(datatable, sep='\s+', header=0)
    vsini = data.vsini
    vsini_err = data.vsini_err

    mean_vsini = np.mean(vsini)
    SEM_vsini = np.std(vsini, ddof=1) / np.sqrt(np.size(vsini)) # standard error of the mean
    propagated_error = 1/(np.sum(1/vsini_err))


    return mean_vsini, SEM_vsini