try:
    import config
    print("Config imported success")
    print(f"Transformers patch applied? {config.transformers.utils.import_utils.check_torch_load_is_safe}")
except Exception as e:
    print(f"Config import failed: {e}")
