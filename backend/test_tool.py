import os
import google.generativeai as genai
import google.ai.generativelanguage as protos
from dotenv import load_dotenv
import json

load_dotenv()
genai.configure(api_key=os.environ['GOOGLE_API_KEY'], transport="rest")

def my_tool(x: int):
    return x * 2

model = genai.GenerativeModel('gemini-3-flash-preview', tools=[my_tool])
chat = model.start_chat()

resp1 = chat.send_message('call my_tool with 5')
print("Model Response 1:", resp1.text if hasattr(resp1, 'text') else "No Text")

# Bypassing the thought_signature error by simply avoiding protos.FunctionResponse
# and just providing the tool result as plain text!
responses_to_send = f"Tool my_tool executed successfully. Result:\n```json\n{json.dumps({'result': 10})}\n```"

print("\nSending plain text instead of FunctionResponse...")
resp2 = chat.send_message(responses_to_send)
print("Model Response 2:")
print(resp2.text)
